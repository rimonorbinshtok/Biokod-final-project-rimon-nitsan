import re
import csv

File_lung = open("data/cBioPorta_TP53_LungCnacer_table.tsv","r")
File_breast = open("data/cBioPortal_breast_cancer_mutations.tsv","r")
output_file = open("results/output_file.txt","w")

def Get_position(Protein_Change):
  match = re.search(r'\d+', Protein_Change)
  if match:
    return int(match.group())
  else:
    print (Protein_Change)
  return None

def count_position(rows, sample_index, mutation_index):
  position_count = {}
  counted_positions = set()
  for row in rows:
    sample_id = row[sample_index]
    protein_change = row[mutation_index]
    position = Get_position(protein_change)
    if position is not None:
      key = (sample_id, position)
      if key not in counted_positions:
        counted_positions.add(key)
        if position not in position_count:
          position_count[position] = 0
        position_count[position] += 1
  return position_count

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

def write_table(output_file, lung_position_count, breast_position_count):
  lung_per = per(lung_position_count)
  breast_per = per(breast_position_count)
  prot_len = 393
  for position in range (1,prot_len+1):
    lung_count = lung_position_count.get(position,0)
    breast_count = breast_position_count.get(position,0)
    lung_percent = lung_per.get(position,0)
    breast_percent = breast_per.get(position, 0)
    output_file.write(f"{position}\t" f"{lung_count}\t"f"{lung_percent:.2f}\t"f"{breast_count}\t"f"{breast_percent:.2f}\n")

### main program
reader_lungs = csv.reader(File_lung, delimiter="\t")
rows_lungs = list(reader_lungs)

reader_breast = csv.reader(File_breast, delimiter="\t")
rows_breast = list(reader_breast)
lung_position_count = count_position(rows_lungs,2 ,5)
breast_position_count = count_position(rows_breast, 2, 5)

write_table(output_file, lung_position_count, breast_position_count)
print("table created")

#test_dict = {
 #   273: 267,
  #  248: 174,
   # 175: 109,
    #245: 132
#}
#print(top(test_dict))
''' Domain_start = 102
Domain_end = 292

def Check_position_in_domain(Protein_Change):
  position = Get_position(Protein_Change)
  if position is None:
    return None
  if Domain_start <= position <= Domain_end:
    return position
  return None
  import re
import csv

File_lung = open("data/cBioPorta_TP53_LungCnacer_table.tsv","r")

def Get_position(Protein_Change):
  match = re.match(r'[A-Z](\d+)[A-Z]', Protein_Change)
  if match:
    return int(match.group(1))
  return None


def count_position_and_mutation_common(rows):
  position_count = {}
  mutation_count = {}
  mutation_position_count = {}
  mutation_type_count = {}
  mutation_type_by_mutations_count ={}
  counted_positions = set()
  counted_mutations = set()
  counted_mutation_types = set()
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
      mutation_type = row[10]
      if mutation_key not in counted_mutations:
        counted_mutations.add(mutation_key)
        if protein_change not in mutation_count:
          mutation_count[protein_change] = 0
        mutation_count[protein_change] += 1
        if position not in mutation_position_count:
          mutation_position_count[position] = 0
        mutation_position_count[position] += 1
        if mutation_type not in mutation_type_by_mutations_count:
          mutation_type_by_mutations_count[mutation_type] = 0
        mutation_type_by_mutations_count[mutation_type] += 1
   #   mutation_type = row[10]
      type_key = (sample_id, mutation_type)
      if type not in counted_mutation_types:
        counted_mutation_types.add(type_key)
        if mutation_type not in mutation_type_count:
          mutation_type_count[mutation_type] = 0
        mutation_type_count[mutation_type]+=1
  return position_count, mutation_count, mutation_position_count, mutation_type_count, mutation_type_by_mutations_count

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
  return pers'''