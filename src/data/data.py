input_file_names = ['IMDB.txt']

def process_file(file_name):
    with open(file_name, 'r') as infile:
        lines = infile.readlines()
    
    with open(file_name, 'w') as outfile:
        for line in lines:
            numbers = line.split()
            if len(numbers) >= 2:
                numbers[0] = str(int(numbers[0]) - 1)
                numbers[1] = str(int(numbers[1]) - 1)
            modified_line = ' '.join(numbers)
            outfile.write(modified_line + '\n')

for file_name in input_file_names:
    process_file(file_name)

print("Processing complete.")
