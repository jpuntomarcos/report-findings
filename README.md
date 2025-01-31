# report-findings
Post-processing code to report events from ENIGMA BRCA1/BRCA2 long read project

### How to use
1. Get Code
```
git clone https://github.com/jpuntomarcos/report-findings.git
```

2. **Configure config.yaml** to define input. Only two parameters require full attention:
   - `ulf_count_paths`: to define paths to Ulf's step 2 output.
   - `events_to_ignore_threshold`. Events equal o smaller than this threshold will be ignored except those already defined in the exons file, like 7_dE7p3
  
3. Launch Script
```
cd report-findings
Rscript report.R config.yaml
```


### Output

The script generates an Excel file containing expected mutations and whether this mutations were found in samples. Columns meaning:

 - FoundInSample: `Yes` if the mutation was found in target sample
 - FoundInControl: `Yes` if the mutation was found in any Control sample
 - FoundInMut: `Yes` if the mutation was found in any Mut sample
 - FIS_n: Number of counts supporting the mutation in target sample
 - FIS_pct: Percentage of counts supporting the mutation in target sample
 - FIS_isoform: Isoform ids (1st part) that were used to compute n.
 - FIS_d: (denominator) Total number of counts that spanned the region of interest. Used to calculate the pct.
 - FIS_d_iso: Isoform ids (1st part) that were used to compute d.
 - FIC_max: Maximum percentage among Control samples. (See FIC column)
 - FIM_max: Maximum percentage among Mut samples. (See FIM column)
 - FIC: Despcription of Control samples where the mutation was found. Within the parenthesis: n and pct. Resuls are sorted by pct, decreasing.
 - FIM: Despcription of Mut samples where the mutation was found. Within the parenthesis: n and pct. Resuls are sorted by pct, decreasing.
 - SampleFindings: Additional finding in the target sample. Describes the full events (i.e. "ins_17:43079334-43079399 del_9") that where identified in the target sample **within the range of interest**. Resuls are sorted by pct, decreasing.
