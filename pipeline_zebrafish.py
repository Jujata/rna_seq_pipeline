#! usr/bin/python3.9

import os
import re
import argparse
import subprocess


def parseArguments():
    """
    Function to parse arguments.

    Input: reads arguments from command line
    Output: parsed arguments from command line (including path)

    """

    # Create the argument parser for RNA-seq analysis
    parser = argparse.ArgumentParser(description="Perform RNA-seq analysis: trimming, mapping, and counting")

    # Define arguments to parse
    parser.add_argument("--input1", "-i", required=True, type=str, help="Path to the first input FASTQ file (read 1)")
    parser.add_argument("--input2", "-I", required=False, type=str, help="Path to the second input FASTQ file (read 2) if pair-end")
    parser.add_argument("--reference_genome", "-r", required=True, type=str, help="Path to the reference genome index")
    parser.add_argument("--annotations", "-a", required=True, type=str, help="Path to the GTF/GFF file with annotations")

    # Call for arguments
    args = parser.parse_args()
    return (args)


def create_main_folder_with_id(input_file):
  """
  Function to extract sample ID and create necessary folders for analysis of FASTQ files (.fastq.gz).

  Args:
    input_file: path to fastq.gz
  """

  # Extract ID samples and create a new folder with its name
  sample_id = re.search(r"([^_]+)(?=\.fastq\.gz))", os.path.basename(input_file)).group(1)
  main_folder = sample_id
  os.makedirs(main_folder, exist_ok=True)

  return main_folder

def fastp(input1, input2, main_folder, output_prefix="trimmed_", cut_tail=20, cut_front=20, cut_mean_quality=20, trim_poly_g=True, trim_poly_x=True, detect_adapter_for_pe=True, html_output="out_FastP.html"):
  """
  This function calls the fastp tool for trimming paired-end reads and moves the outputs to a folder named "trimmed" based on the sample ID.

  Args:
      read1 (str): Path to the read 1 FASTQ file.
      read2 (str): Path to the read 2 FASTQ file.
      output_prefix: Prefix for the output filenames (default: "trimmed_").
      cut_tail: Number of bases to trim from the 3' end (default: 20).
      cut_front: Number of bases to trim from the 5' end (default: 20).
      cut_mean_quality: Minimum average quality score for bases to keep (default: 20).
      trim_poly_g: Trim poly-G tails (default: True).
      trim_poly_x: Trim poly-X tails (default: True).
      detect_adapter_for_pe: Detect adapter sequences for paired-end reads (default: True).
      html_output: Path for the FastQC HTML report (default: "out_FastP.html").
  """

  # Construct output filenames based on prefix and sample ID
  output_r1 = f"{output_prefix}{os.path.splitext(os.path.basename(input1))[0]}.fastq.gz"
  output_r2 = f"{output_prefix}{os.path.splitext(os.path.basename(input2))[0]}.fastq.gz"

  # Create trimmed folder 
  trimmed_folder = os.path.join(main_folder, "trimmed")
  os.makedirs(trimmed_folder, exist_ok=True)

  # Update output filenames with trimmed folder path
  output_r1 = os.path.join(trimmed_folder, output_r1)
  output_r2 = os.path.join(trimmed_folder, output_r2)

  # Build the fastp command as a list of arguments
  fastp_cmd = [
      "fastp",
      "-i", input1,
      "-I", input2,
      "-o", output_r1,
      "-O", output_r2,
      "-q", str(cut_mean_quality),  
      "-l", str(cut_tail),        
      "-u", str(cut_front),       
  ]

  # Add optional arguments based on boolean flags
  if trim_poly_g:
    fastp_cmd.append("--trim_poly_g")
  if trim_poly_x:
    fastp_cmd.append("--trim_poly_x")
  if detect_adapter_for_pe:
    fastp_cmd.append("--detect_adapter_for_pe")

  fastp_cmd.append("-h", html_output)
  fastp_cmd.append(html_output)

  # Call fastp using subprocess
  results_trimming = subprocess.run(fastp_cmd)
  if results_trimming.returncode != 0:
    print("Error running fastp")
    return output_r1, output_r2

  print(f"Finished trimming with fastp. Output files: {output_r1}, {output_r2}")

def hisat2(reference_genome, sample_id, main_folder, input1, input2):
    """
    Function to perform RNA-seq alignment using HISAT2.

    Args:
        reference_genome (str): Path to the reference genome index
        trimmed_folder (str): Path to the folder containing trimmed FASTQ files
        sample_id (str): Sample ID extracted from the filename

    This function uses HISAT2 to align paired-end FASTQ reads against the provided reference genome.
    The output SAM filename is based on the first part of the sample ID (without .fastq.gz).
    """
  # Create a folder for mapping
    mapped_folder = os.path.join(main_folder, "mapped")
    os.makedirs(mapped_folder, exist_ok=True)

    # Construct output SAM filename based on sample name

    output_sam = os.path.join(mapped_folder, f"{sample_id}.sam")

    # Build the HISAT2 command as a list of arguments
    hisat2_cmd = [
        "hisat2",
        "-x", reference_genome,
        "-1", input1,
        "-2", input2,
        "-S", output_sam
    ]

    # Call HISAT2 using subprocess
    results_alignment = subprocess.run(hisat2_cmd)
    if results_alignment.returncode != 0:
        print("Error running HISAT2")
        return output_sam

    print(f"Finished HISAT2 alignment. Output SAM file: {output_sam}")

def samtools(sam_file, main_folder):
    """
    Function to convert a SAM file to a sorted BAM file with indexing and duplicate removal.

    Args:
        sam_file (str): Path to the input SAM file
        main_folder (str): Path to the main folder for output files
    """
    # Define the output BAM file path
    bam_file = os.path.join(main_folder, "mapped", f"{os.path.splitext(os.path.basename(sam_file))[0]}.bam")

    # Convert SAM to unsorted BAM (using -bS)
    samtools_view_cmd = ["samtools", "view", "-bS", sam_file]
    with open(bam_file, 'wb') as bam_out:
        subprocess.run(samtools_view_cmd, stdout=bam_out)

    print(f"Converted {sam_file} to unsorted BAM: {bam_file}")

    # Sort BAM file
    sorted_bam_file = f"{bam_file}.sorted"
    samtools_sort_cmd = ["samtools", "sort", bam_file, "-o", sorted_bam_file]
    subprocess.run(samtools_sort_cmd)
    print(f"Sorted BAM file: {sorted_bam_file}")

    # Index sorted BAM file
    samtools_index_cmd = ["samtools", "index", sorted_bam_file]
    subprocess.run(samtools_index_cmd)
    print(f"Indexed sorted BAM file: {sorted_bam_file}.bai")

    # Remove duplicates (using samtools rmdup)
    dedup_bam_file = os.path.join(main_folder, "mapped", "dedup.bam")
    samtools_rmdup_cmd = ["samtools", "rmdup", sorted_bam_file, dedup_bam_file]
    subprocess.run(samtools_rmdup_cmd)
    print(f"Removed duplicates and created deduplicated BAM: {dedup_bam_file}")

    return dedup_bam_file

def htseq_count(dedup_bam_file, annotation_file, results_folder, main_folder, sample_id):
    """
    Function to run HTSeq-count for gene expression analysis.

    Args:
        dedup_bam_file (str): Path to the deduplicated BAM file
        annotation_file (str): Path to the GTF/GFF annotation file
        results_folder (str): Path to the folder for storing the output file
    """
    # Ensure results folder exists
    results_folder = os.path.join(main_folder, "results")
    os.makedirs(results_folder, exist_ok=True)

    output_results = os.path.join(results_folder, f"{sample_id}.tsv")

    # Construct HTSeq-count command as a list of arguments
    htseq_cmd = [
        "htseq-count",
        "-t", "exon",
        "-i", "gene_id",
        "--stranded=no",
        "-f", "bam",
        "-r", "pos",
        "-s", "no",
        dedup_bam_file,
        annotation_file
    ]

    # Open the output file in write mode
    

    with open(output_results, 'w') as output_file:
        # Run HTSeq-count using subprocess
        results_htseq = subprocess.run(htseq_cmd, stdout=output_file, stderr=subprocess.PIPE)

    if results_htseq.returncode != 0:
        print(f"Error running HTSeq-count: {results_htseq.stderr.decode()}")
        return None  

    print(f"Finished HTSeq-count for {sample_id}. Output: {output_results}")
    return output_results 

if __name__ == "__main__":
    args = parseArguments()
    main_folder = create_main_folder_with_id(args.input1)
    output_r1, output_r2 = fastp(args.input1, args.input2, main_folder)
    if output_r1 and output_r2:
        sam_file = hisat2(args.reference_genome, main_folder, output_r1, output_r2)
        if sam_file:
            dedup_bam_file = samtools(sam_file, main_folder)
            if dedup_bam_file:
                results_folder = os.path.join(main_folder, "results")
                htseq_count(dedup_bam_file, args.annotations, results_folder, os.path.basename(main_folder))
