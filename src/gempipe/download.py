import os
import subprocess


def get_genomes(taxids):
    
    # create a sub-directory
    os.makedirs('working/genomes/', exist_ok=True)
    
    
    # execute the comamnd
    with open('working/genomes/stdout.txt', 'w') as stdout, open('working/genomes/stderr.txt', 'w') as stderr: 
        command = f"""ncbi-genome-download \
            --no-cache \
            --metadata-table working/genomes/metadata.txt \
            --retries 100 --parallel 10 \
            --output-folder working/genomes/ \
            --species-taxids {taxids} \
            --formats assembly-stats,fasta \
            --section genbank \
            bacteria"""
        process = subprocess.Popen(command, shell=True, stdout=stdout, stderr=stderr)
        process.wait()
    
    