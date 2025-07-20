import os
import subprocess
import glob
import csv
import signal
import subprocess
# Base paths
data_base_path = '/data/bipartite_graphs'  # Adjust this to your actual base data path
output_base_path = 'experiments'  # Directory to store outputs
os.makedirs(output_base_path, exist_ok=True)

# Command template
command_template = '../bin/run -f "{data_file}" -p {p} -q {q} -noref'


# Timeout in seconds (e.g., 4 hours)
timeout_seconds = 24 * 60 * 60

# p and q combinations
combinations = [(5, 6), (5, 8), (5, 10), (7, 6), (7, 10), (9, 6), (9, 10)]

# Number of runs per combination
num_runs = 3

# Summary CSV file
summary_csv = 'summary_Norefinement.csv'

# Specific file names to process
specific_file_names = [
    "rec-amz-Books.txt",
    "rec-amz-CDs-and-Vinyl.txt",
    "aff-github-user2project.txt",
    "rec-amz-Apps_for_Android.txt",
    "ia-stackexch-user-marks-post-und.txt",
    "rec-amz-Health-Personal-Care.txt",
    "rec-amz-Grocery-Gourmet-Food.txt"
]

# Function to run a command with timeout and ensure the process group is killed
def run_command(cmd, timeout, output_file):
    timed_out = False
    with open(output_file, "w") as out:
        proc = subprocess.Popen(
            cmd,
            shell=True,
            stdout=out,
            stderr=subprocess.STDOUT,
            start_new_session=True  # Creates a new process group
        )

    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            proc.wait(10)
        except Exception:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)

        with open(output_file, "a") as out:
            out.write("timed out\n")

    return timed_out

# Function to parse the output file
def parse_output(output_content):
    data = {
        'Sample Variance': '',
        'Standard Deviation': '',
        'Biclique Count Final': '',
        'TotalZinShadow Final': '',
        'Time': ''
    }
    lines = output_content.splitlines()
    for line in lines:
        if 'sample variance' in line.lower():
            data['Sample Variance'] = line.split(':')[-1].strip()
        elif 'standard deviation' in line.lower():
            data['Standard Deviation'] = line.split(':')[-1].strip()
        elif 'biclique count final' in line.lower():
            data['Biclique Count Final'] = line.split(':')[-1].strip()
        elif 'totalzinshadow final' in line.lower():
            data['TotalZinShadow Final'] = line.split(':')[-1].strip()
        elif 'time:' in line.lower():
            data['Time'] = line.split(':')[-1].strip()
        elif 'timed out' in line.lower():
            data['Time'] = 'timed out'
    return data

# Main processing loop
results = []

# Create CSV file and write header
with open(summary_csv, 'w', newline='') as csvfile:
    fieldnames = ['Datafile', 'p', 'q', 'Average Count', 'Variance', 'Average Time']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

# Find all .txt files in the base directory and its subdirectories
all_txt_files = glob.glob(os.path.join(data_base_path, '**', '*.txt'), recursive=True)

# Filter files based on the specified file names
files_to_process = [file for file in all_txt_files if os.path.basename(file) in specific_file_names]

if not files_to_process:
    print("No matching files found in the specified directory.")
else:
    for data_file in files_to_process:
        print(f"Processing data file: {data_file}")

        for p, q in combinations:
            counts = []
            times = []
            timeout_occurred = False

            for iteration in range(1, num_runs + 1):
                if timeout_occurred:
                    print(f"Skipping remaining iterations for {data_file}, p={p}, q={q} due to timeout.")
                    break

                output_dir = os.path.join(output_base_path, os.path.basename(data_file), f'p_{p}_q_{q}', f'iteration_{iteration}')
                os.makedirs(output_dir, exist_ok=True)
                output_file = os.path.join(output_dir, 'out.txt')

                cmd = command_template.format(data_file=data_file, p=p, q=q)
                print(f"Running command: {cmd}")
                timeout_occurred = run_command(cmd, timeout_seconds, output_file)

                # Read the output file content
                if os.path.exists(output_file):
                    with open(output_file, 'r') as f:
                        output_content = f.read()
                else:
                    output_content = "timed out"

                # Parse the output
                parsed_data = parse_output(output_content)
                parsed_data['Datafile'] = os.path.basename(data_file)
                parsed_data['p'] = p
                parsed_data['q'] = q

                if 'timed out' in output_content.lower() or timeout_occurred:
                    counts = [-1]
                    times = [-1]
                    break
                else:
                    try:
                        counts.append(float(parsed_data.get('Biclique Count Final', 0)))
                    except ValueError:
                        counts.append(0)
                    try:
                        times.append(float(parsed_data.get('Time', 0)))
                    except ValueError:
                        times.append(0)

                results.append(parsed_data)

            # Calculate average and variance
            if counts and counts[0] != -1:
                avg_count = sum(counts) / len(counts)
                if len(counts) > 1:
                    var_count = sum((x - avg_count) ** 2 for x in counts) / (len(counts) - 1)
                else:
                    var_count = 0.0
                avg_time = sum(times) / len(times)
            else:
                avg_count = var_count = avg_time = -1

            # Write to summary CSV
            with open(summary_csv, 'a', newline='') as csvfile:
                fieldnames = ['Datafile', 'p', 'q', 'Average Count', 'Variance', 'Average Time']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writerow({
                    'Datafile': os.path.basename(data_file),
                    'p': p,
                    'q': q,
                    'Average Count': avg_count,
                    'Variance': var_count,
                    'Average Time': avg_time
                })

print("Processing complete. Summary saved to", summary_csv)