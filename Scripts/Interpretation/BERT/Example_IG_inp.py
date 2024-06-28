# To be used with rxnfp from Scwhaller et al.
# https://github.com/rxn4chemistry/rxnfp.git

from captum.attr import LayerIntegratedGradients
import pandas as pd
from transformers import BertTokenizer, BertModel, BertForMaskedLM, AutoModel, AutoConfig, AutoTokenizer, \
    BertForSequenceClassification, BertPreTrainedModel
import torch.nn as nn
from dotenv import load_dotenv, find_dotenv
import os
import torch
import re

load_dotenv(find_dotenv())
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

    def forward(self, input_ids, attention_mask, token_type_ids, position_ids=None):
        outputs = self.bert(input_ids, attention_mask=attention_mask, token_type_ids=token_type_ids)
        pooled_output = outputs[1]
        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)
        return logits


def get_input_ids(question, tokenizer, device):
    question = SmilesSplit(question)
    question_ids = tokenizer.encode(question, add_special_tokens=True)
    question_ids_without_special_tokens = tokenizer.encode(question, add_special_tokens=False)

    num_tokens = len(question_ids_without_special_tokens)
    reference = '.' * num_tokens
    reference = SmilesSplit(reference)
    reference_ids = tokenizer.encode(reference, add_special_tokens=True)

    question_ids = torch.tensor([question_ids], device=device)
    reference_ids = torch.tensor([reference_ids], device=device)
    return question_ids, reference_ids


def RunIG(question, model, tokenizer):
    question_ids, reference_ids = get_input_ids(question, tokenizer, device)
    seq_length = question_ids.size(1)
    token_type_ids = torch.tensor([[0 for i in range(seq_length)]], device=device)
    attention_mask = torch.ones_like(question_ids)
    question_id_list = question_ids[0].detach().tolist()
    tokens = tokenizer.convert_ids_to_tokens(question_id_list)

    def forward_func(inputs, token_type_ids=None, attention_mask=None):
        pred = model(inputs, token_type_ids=token_type_ids, attention_mask=attention_mask)
        return pred

    question_pred = forward_func(question_ids, token_type_ids=token_type_ids)
    baseline_pred = forward_func(reference_ids, token_type_ids= token_type_ids)

    lig = LayerIntegratedGradients(forward_func, model.bert.embeddings)
    attributions, delta = lig.attribute(inputs=question_ids, baselines=reference_ids, additional_forward_args=(
        token_type_ids, attention_mask), return_convergence_delta=True, n_steps=500)
    attributions_sum = attributions.sum(dim=-1).squeeze(0)

    return tokens, attributions_sum, question_pred.item(), baseline_pred.item()


output_dir = 'BERT_best_hps/BERT_fold_0'
tokenizer_path = 'rxnfp-master/rxnfp/models/transformers/bert_ft/'

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

tokenizer = BertTokenizer.from_pretrained(tokenizer_path, do_basic_tokenize=False, do_lower_case=False)
tokenizer.add_tokens(['[RecipTemp]', '[IonStr]', '[Solv1R]', '[Solv2R]', '0'])

device = torch.device('cuda')
config = AutoConfig.from_pretrained(model_path, num_labels=1)

model = CustomBERTModel.from_pretrained(model_path, config=config)
model.resize_token_embeddings(len(tokenizer))

model.to(device)
model.eval()
model.zero_grad()

full_data = pd.read_excel(r'BERT_input_data.xlsx', sheet_name=None, engine='openpyxl')
test_fold = full_data['Total test fold 0']

with pd.ExcelWriter('IGs_fold_0.xlsx') as writer:
    for i, row in test_fold.iterrows():
           question = row['text']
           tokens, attributions, predictions, baseline_predictions = RunIG(question, model, tokenizer)
           attributions = attributions.detach().cpu().numpy()
           df = pd.DataFrame({'Token': tokens, 'IG': attributions, 'Predictions': predictions,
                           'Baseline predictions': baseline_predictions})
           df.to_excel(writer, sheet_name=f'Total test index {i}', engine='openpyxl')
