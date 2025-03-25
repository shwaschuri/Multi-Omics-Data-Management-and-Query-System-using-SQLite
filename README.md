Overview

This Python script manages and queries multi-omics datasets (metabolomics, proteomics, and transcriptomics) using an SQLite database. It allows users to create the database schema, load multi-omics data from TSV and CSV files, and perform specific queries to analyze subjects and biological entities.
Features

    Database Creation: Automatically creates tables for Metabolomics, Proteomics, Transcriptomics, and Subject metadata.

    Data Loading: Parses and inserts multi-omics data from input files into the SQLite database.

    Querying: Provides predefined queries to extract insights, such as:

        Subjects older than 70 years

        Female subjects with normal BMI

        Visit IDs for a specific subject

        Insulin-resistant subjects with metabolomics data

        Pathways with at least 10 annotations

        Maximum transcript abundance for a specific gene and subject

        Visualization of Age vs. BMI

Requirements

Ensure you have the following dependencies installed before running the script:

    Python 3

    SQLite3

    Pandas

    Matplotlib

Install missing dependencies using:

pip install pandas matplotlib

Usage

Run the script using the following command:

python final.py <command> <database_file>

Commands:
Command	Description
--createdb	Creates the SQLite database with required tables.
--loaddb	Loads metabolomics, proteomics, transcriptomics, and subject metadata into the database.
--querydb=1	Retrieves subjects older than 70 years.
--querydb=2	Finds female subjects with a BMI in the normal range.
--querydb=3	Lists Visit IDs for a given subject.
--querydb=4	Identifies insulin-resistant subjects with metabolomics samples.
--querydb=5	Fetches KEGG IDs for specified metabolite peaks.
--querydb=6	Calculates min, max, and average age of subjects.
--querydb=7	Lists pathways with at least 10 metabolite annotations.
--querydb=8	Finds the maximum abundance of a transcript for a subject.
--querydb=9	Generates a scatter plot of Age vs. BMI.

Example usage:

python final.py --createdb my_database.db
python final.py --loaddb my_database.db
python final.py --querydb=1 my_database.db

Notes

    Ensure the input files (HMP_metabolome_abundance.tsv, HMP_proteome_abundance.tsv, etc.) are present in the working directory.

    Some queries require valid subject IDs and entity IDs present in the database.

    Use python final.py Y my_database.db to check existing tables in the database.

License

This project is open-source and free to use.
