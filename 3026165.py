import re  # import regular expressions module for pattern matching
import sys  # import system-specific parameters and functions
import sqlite3  # import SQLite database module
import pandas as pd  # import pandas for data manipulation and analysis
import matplotlib.pyplot as plt  # import matplotlib for plotting

# Function to parse the sample ID into Subject and Visit components
def sample_parse(sample):
    pattern = r'^([A-Z0-9]+)-([A-Za-z0-9_]+)$'  # regex pattern for parsing sample ID
    match = re.match(pattern, sample)  # match the sample ID to the pattern
    if match:  # if a match is found
        Subject = match.group(1)  # extract subject part
        Visit = match.group(2)  # extract visit part
        return Subject, Visit  # return both components
    return None, None  # return None if no match is found

# Check if the correct number of arguments are provided
if len(sys.argv) < 3:  # ensure the user provides enough arguments
    print("Usage: python final.py <command> <database_file>")  # display usage info
    sys.exit(1)  # exit the script

# Function to clean the metabolite name (removes numbers in parentheses)
def clean_metabolite_name(name):
    cleaned_name = re.sub(r'\(\d+\)$', '', name)  # remove numbers in parentheses at the end
    return cleaned_name.strip()  # return cleaned name without extra spaces

# Connect to the SQLite database
database_file = sys.argv[2]  # get database file from arguments
try: #error handling
    connection = sqlite3.connect(database_file)  # connect to the SQLite database
    print("Database connection successful.")
except sqlite3.Error as e:
    print(f"Error connecting to database: {e}")
    sys.exit(1)
cur = connection.cursor()  # create a cursor to execute SQL commands
if sys.argv[1]=='Y':# to check if the tables are present
    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cur.fetchall()
    print("Tables in the database:", tables)

# Create the necessary tables if the '--createdb' argument is passed
if sys.argv[1] == '--createdb':
    # Create tables for different datasets if they don't already exist
    cur.execute('''
    CREATE TABLE IF NOT EXISTS Metabolite (
        "PeakID" TEXT,
        "Name" TEXT,
        "KEGG_ID" TEXT,
        "HMDB_ID" TEXT,
        "Class" TEXT,
        "Pathway" TEXT
    )
    ''')  # table for metabolite information

    cur.execute('''
    CREATE TABLE IF NOT EXISTS Metabolomics (
        SampleID TEXT NOT NULL,
        PeakID TEXT NOT NULL,
        Abundance REAL,
        FOREIGN KEY (SampleID) REFERENCES Sample(SampleID)
    )
    ''')  # table for metabolomics data

    cur.execute('''
    CREATE TABLE IF NOT EXISTS PeakMetaboliteLink (
        PeakID TEXT NOT NULL,
        MetaboliteID TEXT NOT NULL,
        FOREIGN KEY (PeakID) REFERENCES Metabolomics(PeakID),
        FOREIGN KEY (MetaboliteID) REFERENCES Metabolite(Name)
    )
    ''')  # table linking peaks to metabolites

    cur.execute('''
    CREATE TABLE IF NOT EXISTS Proteomics (
        "SampleID" TEXT NOT NULL,
        "EntityID" TEXT NOT NULL,
        "Abundance" REAL,
        "VisitID" INTEGER,
        FOREIGN KEY("SampleID") REFERENCES "Sample"("SampleID")
    )
    ''')  # table for proteomics data

    cur.execute('''
    CREATE TABLE IF NOT EXISTS Sample (
        "SampleID" TEXT,
        "SubjectID" TEXT NOT NULL,
        "VisitID" TEXT NOT NULL,
        FOREIGN KEY("SubjectID") REFERENCES "Subject"("SubjectID")
    )
    ''')  # table for samples and their relationships to subjects and visits

    cur.execute('''
    CREATE TABLE IF NOT EXISTS "Subject" (
        "Subject_ID" TEXT,
        "Race" TEXT,
        "Sex" TEXT,
        "AGE" NUMERIC,
        "BMI" REAL,
        "SSPG" REAL,
        "IRIS" TEXT
    )
    ''')  # table for subject metadata

    cur.execute('''
    CREATE TABLE IF NOT EXISTS Transcriptomics (
        SampleID TEXT NOT NULL,
        EntityID TEXT NOT NULL,
        Abundance REAL,
        FOREIGN KEY (SampleID) REFERENCES Sample(SampleID)
    )
    ''')  # table for transcriptomics data
    connection.commit()  # save changes to the database

# Load data into the database if the '--loaddb' argument is passed
elif sys.argv[1] == '--loaddb':
    # open and load the metabolomics data
    with open('HMP_metabolome_abundance.tsv', 'r') as meta:
        # read the header of the metabolomics file and split by tabs
        meta_header = meta.readline().strip().split('\t')

        # loop through each line in the metabolomics file
        for line in meta:
            # split the current line into columns by tab
            metabundance = line.strip().split('\t')
            Sample_ID = metabundance[0]  # first column is Sample ID
            SubjID, VisID = sample_parse(Sample_ID)  # parse the sample ID into subject and visit

            # skip invalid sample IDs
            if SubjID is None or VisID is None:
                print(f"Skipping invalid SampleID: {Sample_ID}")
                continue

            # loop through the metabolomics data (from second column onward)
            for i in range(1, len(metabundance)):
                Peak_ID = meta_header[i]  # get the peak ID from header
                Abundance = metabundance[i]  # get the abundance value

                # insert the metabolomics data into the Metabolomics table
                cur.execute('''
                    INSERT INTO Metabolomics (SampleID, PeakID, Abundance)
                    VALUES (?, ?, ?)
                ''', (Sample_ID, Peak_ID, Abundance))

                # insert sample data into the Sample table, avoiding duplicates
                cur.execute('''
                    INSERT OR IGNORE INTO Sample (SampleID, SubjectID, VisitID)
                    VALUES (?, ?, ?)
                ''', (Sample_ID, SubjID, VisID))

    # open and load the subject metadata
    with open('Subject.csv', 'r') as sub:
        header = sub.readline()  # read the header of the subject file
        for line in sub:
            line = line.strip()  # clean the line
            if line:
                try:
                    # split the line by comma and extract values
                    subID, race, sex, age, BMI, sspg, Iris = line.split(',')

                    # insert subject data into the Subject table
                    sql = """
                    INSERT INTO Subject (Subject_ID, Race, Sex, AGE, BMI, SSPG, IRIS)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """
                    cur.execute(sql, (subID, race, sex, age, BMI, sspg, Iris))
                except ValueError:
                    # handle the case of invalid lines
                    print(f"Invalid line: {line}")

    # open and load the proteomics data
    with open('HMP_proteome_abundance.tsv', 'r') as prot:
        # read the header of the proteomics file
        prot_header = prot.readline().strip().split('\t')
        for line1 in prot:
            # split the current line into columns by tab
            protabundance = line1.strip().split('\t')
            Sample_ID1 = protabundance[0]  # first column is Sample ID
            SubjID1, VisID1 = sample_parse(Sample_ID1)  # parse the sample ID into subject and visit

            # skip invalid sample IDs
            if SubjID1 is None or VisID1 is None:
                print(f"Skipping invalid SampleID: {Sample_ID1}")
                continue

            # loop through the proteomics data (from second column onward)
            for i in range(1, len(protabundance)):
                Entity_ID1 = prot_header[i]  # get the entity ID from header
                Abundance1 = protabundance[i]  # get the abundance value

                # insert proteomics data into the Proteomics table
                cur.execute('''
                            INSERT INTO Proteomics (SampleID, EntityID, Abundance, VisitID)
                            VALUES (?, ?, ?, ?)
                            ''', (Sample_ID1, Entity_ID1, Abundance1, VisID1))

                # insert sample data into the Sample table, avoiding duplicates
                cur.execute('''
                    INSERT OR IGNORE INTO Sample (SampleID, SubjectID, VisitID)
                    VALUES (?, ?, ?)
                ''', (Sample_ID1, SubjID1, VisID1))

    # open and load the transcriptomics data
    with open('HMP_transcriptome_abundance.tsv', 'r') as tran:
        # read the header of the transcriptomics file
        tran_header = tran.readline().strip().split('\t')
        for line2 in tran:
            # split the current line into columns by tab
            tranabundance = line2.strip().split('\t')
            Sample_ID2 = tranabundance[0]  # first column is Sample ID
            SubjID2, VisitID2 = sample_parse(Sample_ID2)  # parse the sample ID into subject and visit

            # skip invalid sample IDs
            if SubjID2 is None or VisitID2 is None:
                print(f"Skipping invalid SampleID: {Sample_ID2}")
                continue

            # loop through the transcriptomics data (from second column onward)
            for i in range(1, len(tranabundance)):
                Entity_ID2 = tran_header[i]  # get the entity ID from header
                Abundance2 = tranabundance[i]  # get the abundance value

                # insert transcriptomics data into the Transcriptomics table
                cur.execute('''
                            INSERT INTO Transcriptomics (SampleID, EntityID, Abundance)
                            VALUES (?, ?, ?)
                            ''', (Sample_ID2, Entity_ID2, Abundance2))

                # insert sample data into the Sample table, avoiding duplicates
                cur.execute('''
                               INSERT OR IGNORE INTO Sample (SampleID, SubjectID, VisitID)
                               VALUES (?, ?, ?)
                            ''', (Sample_ID2, SubjID2, VisitID2))

    # open and load the metabolomics annotation data
    with open('HMP_metabolome_annotation.csv', 'r') as meta_ann:
        # read the header of the annotation file
        ann_header = meta_ann.readline().strip().split(',')
        for line3 in meta_ann:
            # split the current line into columns by comma
            PeakID, Metabolite, KEGG, HMDB, ChemicalClass, Pathway = line3.strip().split(',')

            # handle multiple metabolites, KEGG IDs, and HMDB IDs by splitting by pipe character
            Metabolites = Metabolite.split('|') if Metabolite else [None]
            KEGG_IDs = KEGG.split('|') if KEGG else [None]
            HMDB_IDs = HMDB.split('|') if HMDB else [None]

            # loop through each metabolite and insert data
            for i, metabolite in enumerate(Metabolites):
                cleaned_metabolite = clean_metabolite_name(metabolite)  # clean metabolite name
                kegg_id = KEGG_IDs[i] if i < len(KEGG_IDs) else None
                hmdb_id = HMDB_IDs[i] if i < len(HMDB_IDs) else None

                # insert metabolite annotation data into the Metabolite table, avoiding duplicates
                cur.execute('''
                    INSERT OR IGNORE INTO Metabolite (PeakID, Name, KEGG_ID, HMDB_ID, Class, Pathway)
                    VALUES (?, ?, ?, ?, ?, ?)
                ''', (PeakID.strip(), cleaned_metabolite, kegg_id, hmdb_id, ChemicalClass, Pathway))

                # insert the linkage between peak and metabolite
                if PeakID and cleaned_metabolite:
                    cur.execute('''
                        INSERT INTO PeakMetaboliteLink (PeakID, MetaboliteID)
                        VALUES (?, ?)
                    ''', (PeakID.strip(), cleaned_metabolite))

    # commit all changes to the database
    connection.commit()


# Execute various queries based on user input
elif sys.argv[1] == '--querydb=1':
    # Query subjects older than 70 years
    cur.execute('''
        SELECT Subject_ID, AGE
        FROM Subject
        WHERE AGE > 70 AND AGE != 'NA'
    ''')
    rows = cur.fetchall()
    if rows:
        print("Subjects with AGE > 70:")
        for row in rows:
            subject_id, age = row
            print(f"Subject_ID: {subject_id}, AGE: {age}")

elif sys.argv[1] == '--querydb=2':
    # Query female subjects with normal BMI
    cur.execute('''
        SELECT Subject_ID
        FROM Subject
        WHERE Sex = 'F' AND BMI BETWEEN 18.5 AND 24.9
        ORDER BY Subject_ID DESC
    ''')
    rows = cur.fetchall()
    if rows:
        print("Subjects with Sex = 'F' and BMI between 18.5 and 24.9:")
        for row in rows:
            subject_id = row[0]
            print(f"Subject_ID: {subject_id}")


elif sys.argv[1] == '--querydb=3':
    cur.execute('''
        SELECT VisitID
        FROM Sample
        WHERE SubjectID = 'ZNQOVZV'
    ''')
    rows = cur.fetchall()

    if not rows:
        print("No Visit IDs found for Subject 'ZNQOVZV'.")
    else:
        print("Visit IDs for Subject 'ZNQOVZV':")
        unique_visit_ids = set()
        for row in rows:
            visit_id = row[0]
            if visit_id not in unique_visit_ids:
                unique_visit_ids.add(visit_id)
                print(visit_id)

elif sys.argv[1] == '--querydb=4':
    cur.execute('''
        SELECT DISTINCT S.SubjectID
        FROM Sample S
        JOIN Subject SU ON S.SubjectID = SU.Subject_ID
        JOIN Metabolomics M ON S.SampleID = M.SampleID
        WHERE SU.IRIS = 'IR'
    ''')
    rows = cur.fetchall()

    if not rows:
        print("No insulin-resistant subjects with metabolomics samples found.")
    else:
        print("Distinct insulin-resistant SubjectIDs with metabolomics samples:")
        for row in rows:
            print(row[0])

elif sys.argv[1] == '--querydb=5':
    peaks = [
        'nHILIC_121.0505_3.5',
        'nHILIC_130.0872_6.3',
        'nHILIC_133.0506_2.3',
        'nHILIC_133.0506_4.4'
    ]

    cur.execute('''
        SELECT DISTINCT M.KEGG_ID
        FROM PeakMetaboliteLink PML
        JOIN Metabolite M ON PML.MetaboliteID = M.Name
        WHERE PML.PeakID IN (?, ?, ?, ?)
    ''', peaks)

    rows = cur.fetchall()

    if not rows:
        print("No KEGG IDs found for the specified peaks.")
    else:
        print("Unique KEGG IDs for the specified peaks:")
        for row in rows:
            print(row[0])

elif sys.argv[1] == '--querydb=6':
    cur.execute('''
        SELECT MIN(AGE), MAX(AGE), AVG(AGE)
        FROM Subject
        WHERE AGE NOT LIKE 'NA' AND AGE IS NOT NULL AND AGE != ''
    ''')
    row = cur.fetchone()

    if row:
        min_age, max_age, avg_age = row
        print(f"Minimum Age: {min_age}")
        print(f"Maximum Age: {max_age}")
        print(f"Average Age: {avg_age}")

elif sys.argv[1] == '--querydb=7':
    cur.execute('''
        SELECT Pathway, COUNT(*) AS AnnotationCount
        FROM Metabolite
        WHERE Pathway IS NOT NULL AND Pathway != ''
        GROUP BY Pathway
        HAVING COUNT(*) >= 10
        ORDER BY AnnotationCount DESC
    ''')

    rows = cur.fetchall()

    if rows:
        print("Pathways with at least 10 annotations:")
        for row in rows:
            pathway, count = row
            print(f"{pathway}: {count} annotations")

elif sys.argv[1] == '--querydb=8':
    cur.execute('''
        SELECT MAX(T.Abundance) AS MaxAbundance
        FROM Transcriptomics T
        JOIN Sample S ON T.SampleID = S.SampleID
        WHERE S.SubjectID = ? AND T.EntityID = ?;
    ''', ('ZOZOW1T', 'A1BG'))

    result = cur.fetchone()

    if result and result[0] is not None:
        print(f"The maximum abundance of transcript 'A1BG' for subject 'ZOZOW1T' is: {result[0]}")
    else:
        print("No abundance data found for transcript 'A1BG' for subject 'ZOZOW1T'.")

elif sys.argv[1] == '--querydb=9':
    # Plot Age vs. BMI
    query = '''
    SELECT AGE, BMI
    FROM Subject
    WHERE AGE IS NOT NULL AND BMI IS NOT NULL AND AGE != 'NA';
    '''
    df = pd.read_sql_query(query, connection)  # load query result into a DataFrame
    df['AGE'] = pd.to_numeric(df['AGE'], errors='coerce')  # convert AGE to numeric
    df['BMI'] = pd.to_numeric(df['BMI'], errors='coerce')  # convert BMI to numeric
    df = df.dropna()  # drop rows with missing data
    connection.close()# close the database connection
    ax = df.plot.scatter(x='AGE', y='BMI', figsize=(8, 6), color='blue', alpha=0.7, edgecolor='k')  # create scatter plot
    ax.set_title('Scatter Plot of Age vs BMI', fontsize=14)  # set plot title
    ax.set_xlabel('Age', fontsize=12)  # label x-axis
    ax.set_ylabel('BMI', fontsize=12)  # label y-axis
    plt.savefig('age_bmi_scatterplot.png', dpi=300)  # save plot as an image
    plt.show()  # display plot


