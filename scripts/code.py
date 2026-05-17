import re
import csv

File_lung = open("data/cBioPorta_TP53_LungCnacer_table.tsv","r")
File_breast = open("data/bioportal breast cancer.tsv","r")
breast_top = open("results/breast_top.txt","w")
lung_top = open("results/lung_top.txt","w")
output_file = open("results/output_file.txt","w")

def Get_position(Protein_Change): #function to get the mutation position
  match = re.search(r'\d+', Protein_Change)
  if match:
    return int(match.group())
  else:
    return None


def count_position(rows, sample_index, mutation_index): #function that receives the sample ID and the protein change index, and the data from the file. The function returns a dict with the position and amount of samples with a mutation there.
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

def top(count_dict, n=10): # the function receives a dict and the number(n) The function will return a dict with the n highest values from the original dict
  count_dict = count_dict.copy()
  top = {}
  for i in range(n):
    if len(count_dict)== 0:
      break
    max_key = max(count_dict, key=count_dict.get)
    top [max_key] = count_dict[max_key]
    del count_dict[max_key]
  return top

def per(count_dict): # a function that receives a dict and returns a new one where each value represents the percentage of frequency
  pers = {}
  tot = 0
  for key in count_dict:
    tot  += count_dict[key]
  for key in count_dict:
    pers [key] = (count_dict[key]/tot) * 100
  return pers

def write_table(output_file, lung_position_count, breast_position_count): # the function receives 2 dicts (for each cancer type) and writes into a file for each cancer and position the number of samples with a mutation there and the percentage of frequency 
  lung_per = per(lung_position_count)
  breast_per = per(breast_position_count)
  prot_len = 393
  output_file.write("Position\tPatients lung\tPercentage lung\tPatients breast\tPrencentage breast\n")
  for position in range (1,prot_len+1):
    lung_count = lung_position_count.get(position,0)
    breast_count = breast_position_count.get(position,0)
    lung_percent = lung_per.get(position,0)
    breast_percent = breast_per.get(position, 0)
    output_file.write(f"{position}\t" f"{lung_count}\t"f"{lung_percent:.2f}\t"f"{breast_count}\t"f"{breast_percent:.2f}\n")

def write_top_positions(file, top_positions, position_count):#function that writes the top positions, the amount and their percentages to a file 
  file.write("Positions\tPaitents\tPercent\n")
  for position in top_positions:
    file.write(f"{position}\t" f"{position_count[position]}\t" f"{top_positions[position]:.2f}\n")
### main program
reader_lungs = csv.reader(File_lung, delimiter="\t")
rows_lungs = list(reader_lungs)

reader_breast = csv.reader(File_breast, delimiter="\t")
rows_breast = list(reader_breast)
lung_position_count = count_position(rows_lungs,2 ,5)
breast_position_count = count_position(rows_breast, 2, 5)
per_lung_position_count = per(lung_position_count)
per_breast_position_count = per(breast_position_count)
top_lung = top(per_lung_position_count)
top_breast = top(per_breast_position_count)
write_top_positions(breast_top, top_breast, breast_position_count)
write_top_positions(lung_top, top_lung, lung_position_count)

write_table(output_file, lung_position_count, breast_position_count)
print("table created")

breast_top.close()
lung_top.close()
output_file.close()
File_breast.close()
File_lung.close()
  