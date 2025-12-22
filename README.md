```
mermaid
flowchart TD
    Start([Start Pipeline]) --> Config[Configuration & Setup]
    
    Config --> ExcelCheck{Excel File Provided?}
    ExcelCheck -- Yes --> ExtractSRA[Extract SRA IDs from Excel]
    ExcelCheck -- No --> ManualSRA[Manual SRA List Input]
    
    ExtractSRA --> LoadIDs[Load SRA IDs]
    ManualSRA --> LoadIDs
    
    LoadIDs --> ToolCheck{Tools Installed?}
    ToolCheck -- No --> InstallTools[Install All Tools]
    ToolCheck -- Yes --> RefCheck{References Downloaded?}
    
    InstallTools --> RefCheck
    
    RefCheck -- No --> DownloadRefs[Download Reference Files]
    RefCheck -- Yes --> ProcessLoop[Process Each Sample]
    
    DownloadRefs --> ProcessLoop
    
    ProcessLoop --> SampleCheck{Check Sample Status}
    
    SampleCheck -- Completed --> SkipSample[Skipped]
    SampleCheck -- Partial --> ResumeStep[Resume from Last Step]
    SampleCheck -- Not Started --> StartSample[Start Processing]
    
    SkipSample --> NextSample{More Samples?}
    ResumeStep --> ContinueSteps[Continue from Last Step]
    
    StartSample --> Step1[1. Download SRA File]
    Step1 --> ValidateSRA[Validate SRA File]
    ValidateSRA --> Step2[2. Extract FASTQ Files]
    Step2 --> Step3[3. Run FastQC (Raw)]
    Step3 --> Step4[4. Trim Adapters]
    Step4 --> Step5[5. Run FastQC (Trimmed)]
    Step5 --> Step6[6. Align with STAR]
    Step6 --> Step7[7. Variant Calling]
    
    ContinueSteps --> Step7
    
    subgraph VariantCalling [Variant Calling Strategy]
        DataTypeCheck{RNA-seq or DNA-seq?}
        
        DataTypeCheck -- RNA-seq --> RNAStrategy[RNA-seq: Mutect2]
        RNAStrategy --> CreateIntervals[Create Coverage Intervals<br>≥10x coverage]
        CreateIntervals --> RNAMutect2[Run Mutect2]
        
        DataTypeCheck -- DNA-seq --> DNAStrategy[DNA-seq: HaplotypeCaller]
        DNAStrategy --> HaplotypeCaller[Run HaplotypeCaller with GVCF]
    end
    
    Step7 --> VariantCalling
    RNAMutect2 --> Step8[8. Variant Filtration]
    HaplotypeCaller --> Step8
    
    Step8 --> Step9[9. Extract PASS Variants]
    Step9 --> Step10[10. Annotation]
    
    subgraph Annotation [Annotation Strategy]
        FuncotatorCheck{Funcotator Available?}
        
        FuncotatorCheck -- Yes --> UseFuncotator[Use Funcotator]
        FuncotatorCheck -- No --> SnpEffCheck{SnpEff Available?}
        
        SnpEffCheck -- Yes --> UseSnpEff[Use SnpEff]
        SnpEffCheck -- No --> BCFToolsAnno[Use bcftools csq]
        
        BCFToolsAnno --> SimpleAnno[Simple Gene Annotation]
    end
    
    Step10 --> Annotation
    UseFuncotator --> Step11[11. Add dbSNP IDs]
    UseSnpEff --> Step11
    SimpleAnno --> Step11
    
    Step11 --> Step12[12. Add Gene Annotations]
    Step12 --> Step13[13. Export to TSV]
    Step13 --> Cleanup{Keep Intermediate Files?}
    
    Cleanup -- No --> RemoveIntermed[Remove Intermediate Files]
    Cleanup -- Yes --> KeepFiles[Keep All Files]
    
    RemoveIntermed --> NextSample
    KeepFiles --> NextSample
    
    NextSample -- Yes --> ProcessLoop
    NextSample -- No --> MergeTSV[Merge All TSV Files]
    
    MergeTSV --> QCReport[Generate QC Reports]
    QCReport --> FinalOutput[Final Output: all_samples.tsv.gz]
    FinalOutput --> End([Pipeline Complete])
    
    %% Styling
    classDef process fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    classDef decision fill:#f3e5f5
```
