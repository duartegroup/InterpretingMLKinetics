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
import os

load_dotenv(find_dotenv())


torch.manual_seed(42)

SMI_REGEX_PATTERN = r'(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\||\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])'

def SmilesSplit(text):
    regex_pattern = SMI_REGEX_PATTERN
    regex = re.compile(regex_pattern)
    tokens = [token for token in regex.findall(text)]
    text = str(' '.join(tokens))
    return text

def tokenize_function(examples):
    return tokenizer(examples['text'], max_length=300, padding='max_length')

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

def process_outs(outputs, test_data):
    predictions = outputs.predictions
    label_ids = outputs.label_ids
    metrics = outputs.metrics
    predictions_without_brackets = [item[0] for item in predictions]
    prediction_df = pd.DataFrame({'Prediction': predictions_without_brackets, 'Label ID': label_ids.tolist()})
    metric_df = pd.DataFrame(metrics, index=[0])
    prediction_df.to_excel(rf'BERT_fold_0_eval/fold_0_{test_data}_predictions.xlsx')
    metric_df.to_excel(rf'BERT_fold_0_eval/fold_0_{test_data}_metrics.xlsx')


output_dir = 'BERT_best_hps/BERT_fold_0'

epochs = os.listdir(output_dir)
checkpoints = []
checkpoint_numbers = []
for epoch in epochs:
    if 'checkpoint' in epoch:
        checkpoints.append(epoch)
        checkpoint_numbers.append(float(epoch.split('-')[-1]))
    for number, checkpoint in zip(checkpoint_numbers, checkpoints):
        if number == max(checkpoint_numbers):
            final_epoch = checkpoint

model_path = os.path.join(output_dir, final_epoch)
            
tokenizer_path = 'rxnfp-master/rxnfp/models/transformers/bert_ft/'
tokenizer = BertTokenizer.from_pretrained(tokenizer_path, do_basic_tokenize=False, do_lower_case=False)

tokenizer.add_tokens(['[RecipTemp]', '[IonStr]', '[Solv1R]', '[Solv2R]', '0'])
    

full_data = pd.read_excel(r'BERT_input_data.xlsx', sheet_name=None, engine='openpyxl')
train_fold = full_data['Train fold 0'][['text', 'labels']]
test_fold = full_data['Total test fold 0'][['text', 'labels']]

train_fold['text'] = [SmilesSplit(item) for item in train_fold['text']]
test_fold['text'] = [SmilesSplit(item) for item in test_fold['text']]
   

dataset_train = Dataset.from_pandas(train_fold)
dataset_test = Dataset.from_pandas(test_fold)


tokenized_train = dataset_train.map(tokenize_function, batched=True)
tokenized_test = dataset_test.map(tokenize_function, batched=True)
tokenized_train = tokenized_train.remove_columns(['text'])
tokenized_test = tokenized_test.remove_columns(['text'])
tokenized_train.set_format('torch')
tokenized_test.set_format('torch')

training_args = TrainingArguments(output_dir=f'BERT_fold_0_eval')

config = AutoConfig.from_pretrained(model_path, num_labels=1)
model = CustomBERTModel.from_pretrained(model_path, config=config)
model.resize_token_embeddings(len(tokenizer))

device = torch.device('cuda')
model.to(device)
model.eval()
model.zero_grad()


trainer_train = Trainer(model=model, args=training_args,train_dataset = tokenized_train, eval_dataset=tokenized_train)
trainer_test = Trainer(model=model, args=training_args,train_dataset = tokenized_train, eval_dataset=tokenized_test)

outputs_train = trainer_train.predict(tokenized_train)
outputs_test = trainer_test.predict(tokenized_test)


process_outs(outputs_train, 'train')
process_outs(outputs_test, 'total_test')