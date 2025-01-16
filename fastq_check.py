import os # Interact with the file system
import subprocess # Allow the script to use shell commands
import sys # Handle program exit when errors occur

def run_command(command):
    """Run a shell command and handle errors."""
    try:
        subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error while running command: {command}")
        print(e.stderr.decode())
        sys.exit(1)

def check_gzip_integrity(file_path):
    """Check the integrity of a gzipped FASTQ file."""
    print(f"Checking gzip integrity for: {file_path}")
    command = f"gzip -t {file_path}" # Integrity check
    run_command(command) # If there's an error the next print won't run
    print(f"File {file_path} passed gzip integrity check.")

def validate_fastq(file_path):
    """Validate FASTQ format using seqkit."""
    print(f"Validating FASTQ format for: {file_path}")
    command = f"seqkit stats {file_path}"
    run_command(command)
    print(f"FASTQ format validation completed for {file_path}.")

def run_fastqc(file_path, output_dir):
    """Run FastQC to assess read quality."""
    print(f"Running FastQC for: {file_path}")
    command = f"fastqc {file_path} -o {output_dir}"
    run_command(command)
    print(f"FastQC completed for {file_path}. Results saved in {output_dir}.")

def main():
    # Directory containing FASTQ files
    fastq_dir = "path/to/fastq/files"  # Replace with your FASTQ directory
    output_dir = "fastqc_reports"  # Directory for FastQC output
    os.makedirs(output_dir, exist_ok=True)

    # Process each FASTQ file in the directory
    for file_name in os.listdir(fastq_dir):
        file_path = os.path.join(fastq_dir, file_name)

        # Check if the file is gzipped
        if file_name.endswith(".fastq.gz"):
            print(f"Processing file: {file_path}")
            check_gzip_integrity(file_path)  # Step 1: Check gzip integrity
            validate_fastq(file_path)       # Step 2: Validate FASTQ format
            run_fastqc(file_path, output_dir)  # Step 3: Run FastQC
        else:
            print(f"Skipping non-gzipped file: {file_path}")

    print("All FASTQ files have been processed.")

if __name__ == "__main__":
    main()
