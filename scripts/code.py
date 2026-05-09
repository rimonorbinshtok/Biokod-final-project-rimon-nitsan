import re
import csv

File_lung = open("data/cBioPorta_TP53_LungCnacer_table.tsv","r")

def Get_position(Protein_Change):
  match = re.match(r'[A-Z](\d+)[A-Z]', Protein_Change)
  if match:
    return int(match.group(1))
  return None

Domain_start = 102
Domain_end = 292

def Check_position_in_domain(Protein_Change):
  position = Get_position(Protein_Change)
  if position is None:
    return None
  if Domain_start <= position <= Domain_end:
    return position
  return None

def count_position_and_mutation_common(rows):
  position_count = {}
  mutation_count = {}
  mutation_position_count = {}
  counted_positions = set()
  counted_mutations = set()
  #top_10 = {}
  for row in rows:
    sample_id = row[2]
    protein_change = row[5]
    position = Check_position_in_domain(protein_change)
    if position is not None:
      key = (sample_id, position)
      if key not in counted_positions:
        counted_positions.add(key)
        if position not in position_count:
          position_count[position] = 0
        position_count[position] += 1
      mutation_key = (sample_id, protein_change)
      if mutation_key not in counted_mutations:
        counted_mutations.add(mutation_key)
        if protein_change not in mutation_count:
          mutation_count[protein_change] = 0
        mutation_count[protein_change] += 1
        if position not in mutation_position_count:
          mutation_position_count[position] = 0
        mutation_position_count[position] += 1
  return position_count, mutation_count, mutation_position_count

def top(count_dict, n=10):
  count_dict = count_dict.copy()
  top = {}
  for i in range(n):
    if len(count_dict)== 0:
      break
    max_key = max(count_dict, key=count_dict.get)
    top [max_key] = count_dict[max_key]
    del count_dict[max_key]
  return top

def per(count_dict):
  pers = {}
  tot = 0
  for key in count_dict:
    tot  += count_dict[key]
  for key in count_dict:
    pers [key] = (count_dict[key]/tot) * 100
  return pers

### main program

reader = csv.reader(File_lung, delimiter="\t")
rows = list(reader)
position_count, mutation_count, mutation_position_count = count_position_and_mutation_common(rows)
position_per = per(position_count)
mutation_per = per(mutation_count)
mutation_position_per = per(mutation_position_count)
top_positions = top(position_count)
top_mutations = top(mutation_count)
top_positions_by_mutations = top(mutation_position_per)

print("Top positions by patients:")
print(top_positions)
print()
print("Top mutations:")
print(top_mutations)
print()
print("Top positions by mutations:")
print(top_positions_by_mutations)



#test_dict = {
 #   273: 267,
  #  248: 174,
   # 175: 109,
    #245: 132
#}
#print(top(test_dict))
