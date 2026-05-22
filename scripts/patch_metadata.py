import pandas as pd
import numpy as np
import os
import re

def run():
    print("Patching target and sample info metadata...")
    
    # 1. Load target_aadr_metadata.csv as a lookup source
    aadr_meta_path = "/home/duri_bae/popgen_lab/output/target_aadr_metadata.csv"
    if os.path.exists(aadr_meta_path):
        df_aadr = pd.read_csv(aadr_meta_path)
    else:
        print(f"Warning: {aadr_meta_path} not found!")
        df_aadr = pd.DataFrame()
        
    # Clean column names in df_aadr to make them easier to access
    # We can print columns or map them by index
    # Column indices from df_aadr:
    # 0: Genetic ID...
    # 1: Persistent Genetic ID
    # 2: Individual ID
    # 3: Skeletal code
    # 4: Skeletal element (bone)
    # 5: First publication...
    # 6: Publication abbreviation
    # 7: doi
    # 8: Repository
    # 9: Method for Determining Date
    # 10: Date mean in BP
    # 11: Date standard deviation
    # 12: Full Date
    # 13: Age/Sex
    # 14: Group ID
    # 15: Locality
    # 16: Political Entity
    # 17: Latitude
    # 18: Longitude
    # 19: Pulldown Strategy
    # 20: Suffices
    # 21: Data type
    # 22: No. Libraries
    # 23: Mean coverage
    # ...
    # 31: Molecular Sex
    # 35: Terminal Y Haplogroup
    # 38: mtDNA haplogroup
    # 43: Library type (minus/half/plus...)
    
    # Create lookup dictionaries based on AADR metadata
    aadr_lookup = {}
    if not df_aadr.empty:
        for idx, row in df_aadr.iterrows():
            # Standardize string fields
            ind_id = str(row.iloc[2]).strip()
            skeletal_code = str(row.iloc[3]).strip()
            genetic_id = str(row.iloc[0]).strip()
            
            meta = {
                'Cal': str(row.iloc[12]).strip() if not pd.isna(row.iloc[12]) else (str(row.iloc[10]).strip() + ' BP' if not pd.isna(row.iloc[10]) else 'unknown'),
                'Period': str(row.iloc[14]).strip() if not pd.isna(row.iloc[14]) else 'unknown', # Group ID
                'Lifestyle': 'unknown',
                'Culture': 'unknown',
                'Country': str(row.iloc[16]).strip() if not pd.isna(row.iloc[16]) else 'unknown',
                'Paper': str(row.iloc[6]).strip() if not pd.isna(row.iloc[6]) else 'unknown',
                'City': str(row.iloc[15]).strip() if not pd.isna(row.iloc[15]) else 'unknown',
                'Library': str(row.iloc[21]).strip() if not pd.isna(row.iloc[21]) else 'unknown',
                'Data type': str(row.iloc[21]).strip() if not pd.isna(row.iloc[21]) else 'unknown',
                'UDG': str(row.iloc[43]).strip() if not pd.isna(row.iloc[43]) else 'unknown',
                'latitude': str(row.iloc[17]).strip() if not pd.isna(row.iloc[17]) else 'unknown',
                'longitude': str(row.iloc[18]).strip() if not pd.isna(row.iloc[18]) else 'unknown',
                'Published': str(row.iloc[6]).strip() if not pd.isna(row.iloc[6]) else 'unknown',
                'bone': str(row.iloc[4]).strip() if not pd.isna(row.iloc[4]) else 'unknown',
                'source_file': 'v66.1240K.aadr.PUB.csv',
                'match_status': 'iid_match'
            }
            # We can lookup by individual ID
            if ind_id and ind_id != 'nan':
                aadr_lookup[ind_id] = meta.copy()
            # Also lookup by skeletal code
            if skeletal_code and skeletal_code != 'nan':
                meta_skeletal = meta.copy()
                meta_skeletal['match_status'] = 'skeletal_code_match'
                aadr_lookup[skeletal_code] = meta_skeletal
            # Also lookup by genetic ID
            if genetic_id and genetic_id != 'nan':
                meta_genetic = meta.copy()
                meta_genetic['match_status'] = 'genetic_id_match'
                aadr_lookup[genetic_id] = meta_genetic

    # 2. Add SUI001.WGS specifically since we know it maps to U1
    # We found in AADR: U1 has library SUI001.A0101, dates 3375 BP, Mustang, Nepal (29.04, 83.93)
    sui_meta = {
        'Cal': '1550-1300 BCE (3375 BP)',
        'Period': 'Nepal_LateC',
        'Lifestyle': 'unknown',
        'Culture': 'unknown',
        'Country': 'Nepal',
        'Paper': 'LiuJeongNatComm2022',
        'City': 'Suila (Mustang District)',
        'Library': 'ss.minus',
        'Data type': '1240k',
        'UDG': 'minus',
        'latitude': '29.04',
        'longitude': '83.93',
        'Published': 'LiuJeongNatComm2022',
        'bone': 'tooth',
        'source_file': 'v66.1240K.aadr.PUB.csv',
        'match_status': 'manual_u1_match'
    }
    aadr_lookup['SUI001.WGS'] = sui_meta

    # Mbc individual date map for MBC1-9
    mbc_2024_dates = {
        'MBC1': '4,960–4,827 cal BP (mean: 4894 BP)',
        'MBC2': '4,852–4,625 cal BP (mean: 4739 BP)',
        'MBC4': '5,290–4,884 cal BP (mean: 5087 BP)',
        'MBC5': '5,211–4,855 cal BP (mean: 5033 BP)',
        'MBC6': '5,264–4,860 cal BP (mean: 5062 BP)',
        'MBC7': '5,042–4,869 cal BP (mean: 4956 BP)',
        'MBC8': '5,214–4,862 cal BP (mean: 5038 BP)',
        'MBC9': '4,813–4,449 cal BP (mean: 4631 BP)'
    }

    def patch_row(row):
        iid = str(row['IID']).strip()
        pid = str(row['PID']).strip()
        
        # Check if we have a match in AADR lookup
        if iid in aadr_lookup:
            meta = aadr_lookup[iid]
            for col in meta:
                if col in row.index:
                    row[col] = meta[col]
            return row
            
        # Check if Mbc sample
        if iid.startswith('Mabuco') or iid.startswith('MBC'):
            row['latitude'] = '28.31'
            row['longitude'] = '89.43'
            row['Country'] = 'China'
            row['City'] = 'Kangmar County, Shigatse'
            row['Lifestyle'] = 'unknown'
            row['Culture'] = 'unknown'
            
            # Check if 2024 paper sample (MBC1-9)
            if iid in mbc_2024_dates:
                row['Cal'] = mbc_2024_dates[iid]
                row['Period'] = 'Neolithic'
                row['Paper'] = 'Ran2024NatEcolEvol'
                row['Published'] = 'Ran2024NatEcolEvol'
                row['bone'] = 'petrous (cochlea)'
                row['source_file'] = 'Ran2024NatEcolEvol'
                row['match_status'] = 'manual_lookup'
            else:
                # 2025 iScience sample: Mbc4.4k, Mbc4k, Mbc3.5k based on PID suffix
                row['Paper'] = 'Ran2025iScience'
                row['Published'] = 'Ran2025iScience'
                row['bone'] = 'bone'
                row['source_file'] = 'Ran2025iScience'
                row['match_status'] = 'manual_lookup_pid_suffix'
                
                # Check PID suffix for date
                if '4.4k' in pid:
                    row['Cal'] = '4400 BP'
                    row['Period'] = 'LN'
                elif '4k' in pid:
                    row['Cal'] = '4000 BP'
                    row['Period'] = 'LN/EBA'
                elif '3.5k' in pid:
                    row['Cal'] = '3500 BP'
                    row['Period'] = 'BA'
                else:
                    row['Cal'] = '4400-3500 BP'
                    row['Period'] = 'unknown'
            return row
            
        # Fallback for target samples that might still be empty/unknown but we want to fill
        # e.g., C4777 is not in AADR but is Zongri4.5k
        if iid == 'C4777':
            row['latitude'] = '35.3'
            row['longitude'] = '100.4'
            row['Country'] = 'China'
            row['City'] = 'Zongri (Qinghai)'
            row['Cal'] = '3400-1900 BCE'
            row['Period'] = 'Zongri4.5k'
            row['Paper'] = 'WangFuSciAdv2023'
            row['Published'] = 'WangFuSciAdv2023'
            row['bone'] = 'bone'
            row['source_file'] = 'genotype'
            row['match_status'] = 'manual_lookup_zongri'
            return row
            
        return row

    # Patch target_metadata_260522.csv
    target_path = "/home/duri_bae/popgen_lab/output/target_metadata_260522.csv"
    if os.path.exists(target_path):
        df_target = pd.read_csv(target_path)
        # Apply patch row-by-row
        df_target = df_target.apply(patch_row, axis=1)
        df_target.to_csv(target_path, index=False)
        print(f"Patched target_metadata_260522.csv and saved to {target_path}")
        
    # Patch SampleInfo_new_260522.csv
    sample_info_path = "/home/duri_bae/popgen_lab/output/SampleInfo_new_260522.csv"
    if os.path.exists(sample_info_path):
        df_sample = pd.read_csv(sample_info_path)
        df_sample = df_sample.apply(patch_row, axis=1)
        df_sample.to_csv(sample_info_path, index=False)
        print(f"Patched SampleInfo_new_260522.csv and saved to {sample_info_path}")

if __name__ == '__main__':
    run()
