
# To be used with rxnfp from Scwhaller et al.
# https://github.com/rxn4chemistry/rxnfp.git

import pandas as pd
from transformers import BertTokenizer, BertModel, BertForMaskedLM, AutoModel, AutoConfig, AutoTokenizer, \
    BertForSequenceClassification, BertPreTrainedModel
import torch
import torch.nn as nn
import re
from datasets import Dataset, list_metrics, load_metric
from torch.nn import MSELoss
from transformers import TrainingArguments, Trainer
from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

from ax.service.managed_loop import optimize

torch.manual_seed(42)

SMI_REGEX_PATTERN = r'(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\||\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])'


def SmilesSplit(text):
    regex_pattern = SMI_REGEX_PATTERN
    regex = re.compile(regex_pattern)
    tokens = [token for token in regex.findall(text)]
    text = str(' '.join(tokens))
    return text


class CustomBERTModel(BertPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)
        self.num_labels = 1
        self.bert = BertModel(config)
        self.dropout = nn.Dropout(config.hidden_dropout_prob)
        self.classifier = nn.Linear(256, 1)

    def forward(self, input_ids, attention_mask, token_type_ids, labels, position_ids=None):
        outputs = self.bert(input_ids, attention_mask=attention_mask, token_type_ids=token_type_ids)
        pooled_output = outputs[1]
        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)
        loss_fct = MSELoss()
        loss = loss_fct(logits.squeeze(), labels.squeeze())
        return (loss,) + (logits,)


def train_evaluate(parameterization):
    model_path = 'rxnfp-master/rxnfp/models/transformers/bert_ft/'
    tokenizer = BertTokenizer.from_pretrained(model_path, do_basic_tokenize=False, do_lower_case=False)
    tokenizer.add_tokens(['[RecipTemp]', '[IonStr]', '[Solv1R]', '[Solv2R]', '0'])

    config = AutoConfig.from_pretrained(model_path, num_labels=1)
    config.hidden_dropout_prob = parameterization['hidden_dropout_prob']

    def tokenize_function(examples):
        return tokenizer(examples['text'], max_length=300, padding='max_length')

    full_data = pd.read_excel(r'BERT_input_data.xlsx', sheet_name=None, engine='openpyxl')
    train_fold = full_data['Train fold 0'][['text', 'labels']]
    val_fold = full_data['Val fold 0'][['text', 'labels']]

    train_fold['text'] = [SmilesSplit(item) for item in train_fold['text']]
    val_fold['text'] = [SmilesSplit(item) for item in val_fold['text']]

    dataset_train = Dataset.from_pandas(train_fold)
    dataset_test = Dataset.from_pandas(val_fold)

    tokenized_train = dataset_train.map(tokenize_function, batched=True)
    tokenized_test = dataset_test.map(tokenize_function, batched=True)

    tokenized_train = tokenized_train.remove_columns(['text'])
    tokenized_test = tokenized_test.remove_columns(['text'])
    tokenized_train.set_format('torch')
    tokenized_test.set_format('torch')

    training_args = TrainingArguments(output_dir=f"BERT_fold_0_LR_{parameterization['learning_rate']}_D_{parameterization['hidden_dropout_prob']}", evaluation_strategy = 'epoch', per_device_train_batch_size = 16, learning_rate = parameterization['learning_rate'], lr_scheduler_type = 'linear', num_train_epochs = 10, save_strategy = 'epoch')

    model = CustomBERTModel.from_pretrained(model_path, config=config)
    model.resize_token_embeddings(len(tokenizer))
    device = torch.device('cuda')
    model.to(device)

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_train,
        eval_dataset=tokenized_test,
    )

    trainer.train()
    return trainer.evaluate(tokenized_test)['eval_loss']


best_parameters, values, experiment, model = optimize(parameters=[
    {'name': 'learning_rate', 'type': 'range', 'bounds': [1e-6, 1e-4]},
    {'name': 'hidden_dropout_prob', 'type': 'range', 'bounds': [0.05, 0.8]}],
    evaluation_function=train_evaluate,
    objective_name='eval_loss', total_trials=20, minimize=True)

means, covs = values

parameter_df = pd.DataFrame({k: [v] for k, v in best_parameters.items()})
values_df = pd.DataFrame({'Objective mean': [means['eval_loss']], 'SEM squared': [covs['eval_loss']['eval_loss']]})
total_df = pd.concat((parameter_df, values_df), axis=1)
total_df.to_excel(r'hp_sweep_fold_0_results.xlsx', engine='openpyxl')
