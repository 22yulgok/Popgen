import pandas as pd
import numpy as np
import os

def clean_date(date_val):
    if pd.isna(date_val):
        return -1
    try:
        # Some dates might be strings or have weird characters. Just try converting to float.
        # But wait, looking at the dataset some might be ".."
        if str(date_val).strip() == "..":
            return -1
        return float(date_val)
    except ValueError:
        return -1

def run():
    print("Loading data...")
    # Load original sample info
    sample_info_path = "/home/duri_bae/popgen_lab/raw/SampleInfo_260522.csv"
    df_sample = pd.read_csv(sample_info_path)
    # Ensure IID is treated as string
    df_sample['IID'] = df_sample['IID'].astype(str).str.strip()
    
    # Load IND file (space/tab separated, columns typically IID, SEX, PID)
    ind_path = "/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.ind"
    df_ind = pd.read_csv(ind_path, sep=r'\s+', header=None, names=["IID", "SEX", "PID"])
    df_ind['IID'] = df_ind['IID'].astype(str).str.strip()
    df_ind['PID'] = df_ind['PID'].astype(str).str.strip()

    # Load AADR file
    aadr_path = "/home/duri_bae/popgen_lab/raw/v66.1240K.aadr.PUB.csv"
    # AADR has many columns, but let's load it and parse the important ones
    df_aadr = pd.read_csv(aadr_path, dtype=str)
    
    # Identify relevant columns in AADR
    aadr_cols = df_aadr.columns.tolist()
    genetic_id_col = aadr_cols[0]  # Genetic ID...
    ind_id_col = aadr_cols[2]      # Individual ID
    pid_col = aadr_cols[14]        # Group ID
    date_col = aadr_cols[10]       # Date mean in BP...
    pub_col = aadr_cols[6]         # Publication abbreviation
    lat_col = aadr_cols[17]        # Latitude
    lon_col = aadr_cols[18]        # Longitude
    
    # Clean up AADR
    df_aadr['Date_Numeric'] = df_aadr[date_col].apply(clean_date)
    
    # 1. Which IIDs are in IND but not in SampleInfo?
    existing_iids = set(df_sample['IID'].tolist())
    ind_iids = set(df_ind['IID'].tolist())
    missing_iids = ind_iids - existing_iids
    
    # 2. Build mapping of PID -> max date across all AADR
    # This helps with the rule: "Exception: If the same pid contains both samples with dates >= 4000 
    # and samples with dates < 4000, include the samples from that pid as well."
    # Basically, does this PID have ANY sample with date >= 4000?
    pid_max_dates = df_aadr.groupby(pid_col)['Date_Numeric'].max().to_dict()
    
    print(f"Total samples in genotype metadata: {len(ind_iids)}")
    print(f"Samples already present in SampleInfo: {len(ind_iids.intersection(existing_iids))} (out of {len(existing_iids)} total in SampleInfo)")
    print(f"Missing samples identified: {len(missing_iids)}")
    
    # Setup tracking stats
    added_samples = []
    matched_iid = 0
    matched_pid_only = 0
    not_found = 0
    included_by_pid_exception = 0
    
    # Get columns of the original sample info
    output_columns = list(df_sample.columns)
    # We might add tracking columns if desired, but user said "add them only if necessary". We'll just stick to original columns and fill them.
    # Actually user suggested: "If additional columns are needed for tracking, add them only if necessary and clearly name them, for example: source_file, match_status, date_filter_status, notes"
    # Let's add them to all rows. Existing rows will have NA or "original".
    
    df_sample['source_file'] = 'SampleInfo_260522.csv'
    df_sample['match_status'] = 'original'
    df_sample['date_filter_status'] = 'N/A'
    df_sample['notes'] = ''
    
    output_columns.extend(['source_file', 'match_status', 'date_filter_status', 'notes'])
    
    for _, row in df_ind.iterrows():
        iid = row['IID']
        if iid not in missing_iids:
            continue
            
        pid = row['PID']
        sex = row['SEX']
        
        # 3. Match in AADR
        # Try exact IID match first
        # We can check genetic_id_col and ind_id_col
        match = df_aadr[(df_aadr[genetic_id_col] == iid) | (df_aadr[ind_id_col] == iid)]
        
        match_type = 'not_found'
        aadr_row = None
        
        if len(match) > 0:
            match_type = 'iid'
            aadr_row = match.iloc[0]
        else:
            # Try PID match
            pid_match = df_aadr[df_aadr[pid_col] == pid]
            if len(pid_match) > 0:
                match_type = 'pid_only'
                # Just take the first one from that group to fill in "Period", "Culture", etc.
                # BUT wait, "Do not infer sample information unless the match is clearly supported by iid or pid."
                # If we use pid_only, we can take some group-level info, but IID is still not found.
                aadr_row = pid_match.iloc[0]
        
        if match_type == 'not_found':
            # Include anyway, write 'unknown'/'not_found'
            not_found += 1
            new_row = {col: 'not_found' for col in output_columns}
            new_row['IID'] = iid
            new_row['SEX'] = sex
            new_row['PID'] = pid
            new_row['source_file'] = 'genotype'
            new_row['match_status'] = 'not_found'
            new_row['date_filter_status'] = 'N/A'
            added_samples.append(new_row)
            continue
            
        # We have a match. Check the date filter.
        # "Include samples whose value in this column is greater than or equal to 4000."
        date_numeric = aadr_row['Date_Numeric']
        
        passes_date = False
        date_status = ''
        if date_numeric >= 4000:
            passes_date = True
            date_status = 'date >= 4000'
        else:
            # Check the PID exception
            aadr_pid = aadr_row[pid_col]
            if not pd.isna(aadr_pid) and aadr_pid in pid_max_dates and pid_max_dates[aadr_pid] >= 4000:
                passes_date = True
                date_status = 'exception: group pid >= 4000'
                included_by_pid_exception += 1
            else:
                passes_date = False
                
        # Wait, the instruction says:
        # "Include samples whose value... >= 4000"
        # "Exception: If the same pid contains both samples with dates >= 4000 and samples with dates < 4000, include the samples from that pid as well."
        # If it doesn't pass the date rule, do we exclude it?
        # Re-read: "For samples whose information cannot be found in /raw/v66.1240K.aadr.PUB.csv, still include the sample in the output and write not_found..."
        # BUT for samples whose info IS found, but they fail the date filter, what do we do? 
        # Usually we'd exclude them, but the prompt says: "Update the sample information table by identifying samples that are present in the genotype metadata but missing from the current sample information file, then retrieve their metadata from the AADR public annotation file... When retrieving information... apply the following date rule..." 
        # It seems the date rule applies to the retrieval part. Let's assume we ONLY include if it passes the date rule, OR if it's not found at all.
        # Wait, if we exclude it, it means it's not added to the output.
        # Let's read carefully: "In other words, include a missing sample if: - its own date mean BP is >= 4000, or - its pid has at least one associated sample in the AADR file with date mean BP >= 4000."
        # "For samples whose information cannot be found in /raw/v66.1240K.aadr.PUB.csv, still include the sample in the output and write not_found or unknown in the unavailable fields."
        # OK, so if it fails the date filter, we SKIP it entirely and do NOT add it.
        
        if not passes_date:
            continue
            
        # Construct the new row
        if match_type == 'iid':
            matched_iid += 1
        else:
            matched_pid_only += 1
            
        new_row = {col: 'unknown' for col in output_columns}
        new_row['IID'] = iid
        new_row['SEX'] = sex
        new_row['PID'] = pid
        
        # Fill in what we can from aadr_row
        # SampleInfo_260522 columns: IID,SEX,PID,Cal,Period,Lifestyle,Culture,Country,Paper,City,Library,Data type,UDG,latitude,longitude,Published,bone
        new_row['latitude'] = aadr_row[lat_col] if not pd.isna(aadr_row[lat_col]) else 'unknown'
        new_row['longitude'] = aadr_row[lon_col] if not pd.isna(aadr_row[lon_col]) else 'unknown'
        new_row['Paper'] = aadr_row[pub_col] if not pd.isna(aadr_row[pub_col]) else 'unknown'
        new_row['bone'] = aadr_row[aadr_cols[4]] if not pd.isna(aadr_cols[4]) and not pd.isna(aadr_row[aadr_cols[4]]) else 'unknown'
        # Country might be political entity (col 16)
        new_row['Country'] = aadr_row[aadr_cols[16]] if not pd.isna(aadr_row[aadr_cols[16]]) else 'unknown'
        # Cal could be full date (col 12)
        new_row['Cal'] = aadr_row[aadr_cols[12]] if not pd.isna(aadr_row[aadr_cols[12]]) else 'unknown'
        
        new_row['source_file'] = 'v66.1240K.aadr.PUB.csv'
        new_row['match_status'] = match_type
        new_row['date_filter_status'] = date_status
        new_row['notes'] = f"Matched PID {aadr_row[pid_col]}" if match_type == 'pid_only' else ''
        
        added_samples.append(new_row)

    print(f"Missing samples added to the output: {len(added_samples)}")
    print(f"Added samples matched by IID: {matched_iid}")
    print(f"Added samples matched only by PID: {matched_pid_only}")
    print(f"Samples marked as not_found or unknown (due to no match): {not_found}")
    print(f"Samples included because of the same-pid date exception: {included_by_pid_exception}")
    
    # Concatenate and save
    df_added = pd.DataFrame(added_samples)
    if len(df_added) > 0:
        # Ensure same order of columns
        df_added = df_added[output_columns]
        df_final = pd.concat([df_sample, df_added], ignore_index=True)
    else:
        df_final = df_sample
        
    out_dir = "/home/duri_bae/popgen_lab/output"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "SampleInfo_new_260522.csv")
    df_final.to_csv(out_path, index=False)
    print(f"Saved output to {out_path}")

if __name__ == '__main__':
    run()
