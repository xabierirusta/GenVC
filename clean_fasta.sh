#!/bin/bash

# Clean the reference fasta file for CNVkit
sed 's/^>\([^chr][^ ]*\)/>chr\1/' "${params.ref}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" > "${params.ref}/Homo_sapiens.GRCh38.dna.primary_assembly_with_chr.fa"
sed 's/^>.*chromosome:GRCh38:\([0-9XY]*\):.*$/>chr\1/' "${params.ref}/Homo_sapiens.GRCh38.dna.primary_assembly_with_chr.fa" > "${params.ref}/cleaned_reference.fa"
