# Overall Pipeline Architecture

```mermaid
graph TD
    Start[Start: MDD2 Pipeline] --> Config{Configuration Check}
    
    Config --> Excel{Excel File Provided?}
    Excel -->|Yes| ExtractSRA[Extract SRA IDs from Excel]
    Excel -->|No| Prompt[Prompt for Excel File]
    Prompt --> ExtractSRA
    
    ExtractSRA --> SRA_List[List of SRA IDs]
    
    SRA_List --> Command{Command Selection}
    
    Command -->|install| InstallFlow[Tool Installation Flow]
    Command -->|setup| SetupFlow[Reference Setup Flow]
    Command -->|run| MainFlow[Main Analysis Flow]
    Command -->|fastq| FastqFlow[FASTQ Processing Flow]
    Command -->|align| AlignFlow[Alignment Flow]
    Command -->|process-all| ProcessFlow[Variant Calling Flow]
    
    InstallFlow --> End1[End: Tools Installed]
    SetupFlow --> End2[End: References Ready]
    MainFlow --> Final[End: Analysis Complete]
    FastqFlow --> End3[End: FASTQ Processed]
    AlignFlow --> End4[End: Alignment Complete]
    ProcessFlow --> End5[End: Variants Called]
```
# Tool Installation Flow
```mermaid
graph TD
    StartInstall[Start: install] --> SetupEnv[Setup Environment]
    
    SetupEnv --> Java[Install Java]
    Java --> SRA[Install SRA Toolkit]
    SRA --> FastQC[Install FastQC]
    FastQC --> MultiQC[Install MultiQC]
    MultiQC --> Trimmomatic[Install Trimmomatic]
    Trimmomatic --> GATK[Install GATK]
    GATK --> BioTools[Install Bioinformatics Tools<br/>samtools/bcftools]
    BioTools --> Parallel[Install Parallel Tools<br/>pigz/parallel]
    Parallel --> STAR[Install STAR Aligner]
    
    STAR --> Verify[Verify All Installations]
    Verify --> CreateEnv[Create Environment Script]
    CreateEnv --> EndInstall[End: Tools Ready]
```
#  Main Analysis Flow (Sequential Processing)

```mermaid
graph TD
    StartRun[Start: run command] --> ForEach[For Each SRA ID]
    
    ForEach --> CheckStatus{Check Sample Status}
    
    CheckStatus -->|not_started| Download[Download SRA]
    CheckStatus -->|downloaded| Extract[Extract FASTQ]
    CheckStatus -->|extracted| Trim[Trim Reads]
    CheckStatus -->|trimmed| Align[Align with STAR]
    CheckStatus -->|aligned| Process[Process Variants]
    CheckStatus -->|completed| Skip[Skip - Already Done]
    
    Download --> QC1[FastQC: Raw Reads]
    Extract --> QC1
    QC1 --> Trim
    Trim --> QC2[FastQC: Trimmed Reads]
    QC2 --> Align
    Align --> Process
    Process --> NextSample{More Samples?}
    
    NextSample -->|Yes| ForEach
    NextSample -->|No| Merge[Merge All TSV Files]
    Merge --> FinalOutput[Final: all_samples.tsv.gz]
    FinalOutput --> EndRun[End: Pipeline Complete]
```

# Single Sample Processing Pipeline

```mermaid
graph TD
    StartSample[Sample: SRRXXXXXX] --> DownloadSRA[Download SRA File]
    DownloadSRA --> ValidateSRA[Validate SRA]
    ValidateSRA --> ExtractFASTQ[Extract Paired FASTQ]
    ExtractFASTQ --> RawQC[FastQC on Raw Reads]
    
    RawQC --> TrimAdapters[Trim with Trimmomatic]
    TrimAdapters --> TrimmedQC[FastQC on Trimmed Reads]
    
    TrimmedQC --> CheckSTAR{STAR Index Exists?}
    CheckSTAR -->|No| BuildIndex[Build STAR Index]
    CheckSTAR -->|Yes| AlignSTAR[Align with STAR]
    BuildIndex --> AlignSTAR
    
    AlignSTAR --> IndexBAM[Index BAM File]
    
    IndexBAM --> ProcessDNA{RNA or DNA?}
    
    ProcessDNA -->|DNA| DNA_Flow[DNA Processing Flow]
    ProcessDNA -->|RNA| RNA_Flow[RNA Processing Flow]
    
    DNA_Flow --> TSV_DNA[Export DNA Variants to TSV]
    RNA_Flow --> TSV_RNA[Export RNA Variants to TSV]
    
    TSV_DNA --> Cleanup{Keep Intermediate?}
    TSV_RNA --> Cleanup
    
    Cleanup -->|Yes| Keep[Keep All Files]
    Cleanup -->|No| Remove[Remove Intermediate Files]
    
    Keep --> EndSample[End: Sample Complete]
    Remove --> EndSample
```
# DNA-seq Variant Calling Flow

```mermaid
graph TD
    StartDNA[Input: Aligned BAM] --> AddRG[Add/Replace Read Groups]
    
    AddRG --> MarkDup[Mark Duplicates]
    MarkDup --> SplitNCigar[Split N Cigar Reads]
    SplitNCigar --> BQSR1[BaseRecalibrator<br/>Create Recal Table]
    BQSR1 --> ApplyBQSR[Apply BQSR]
    
    ApplyBQSR --> CallVariants[HaplotypeCaller<br/>ERC GVCF Mode]
    CallVariants --> Filter[Variant Filtration]
    
    Filter --> Annotate{Annotation Method}
    
    Annotate -->|Funcotator| Func[Funcotator Annotation]
    Annotate -->|SnpEff| SnpEff[SnpEff Annotation]
    Annotate -->|Fallback| Simple[Simple Annotation]
    
    Func --> AddRSID[Add dbSNP IDs]
    SnpEff --> AddRSID
    Simple --> AddRSID
    
    AddRSID --> AddGenes[Add Gene Annotations]
    AddGenes --> Export[Export to TSV]
    Export --> EndDNA[End: DNA Variants]
```
# RNA-seq Variant Calling Flow Evidence-Based
```mermaid
graph TD
    StartRNA[Input: RNA-seq BAM] --> AddRG_RNA[Add Read Groups]
    
    AddRG_RNA --> MarkDup_RNA[Mark Duplicates]
    MarkDup_RNA --> SplitNCigar_RNA[Split N Cigar Reads<br/>RNA-specific]
    
    SplitNCigar_RNA --> BQSR_Opt{BQSR for RNA?}
    
    BQSR_Opt -->|Yes| BQSR_RNA[Base Quality Recalibration]
    BQSR_Opt -->|No| SkipBQSR[Skip BQSR]
    
    BQSR_RNA --> CallVariants_RNA
    SkipBQSR --> CallVariants_RNA
    
    CallVariants_RNA --> Method{RNA Variant Calling Method}
    
    Method -->|Evidence-Based| CoverageBed[Create Coverage BED<br/>=10x coverage]
    CoverageBed --> Mutect2[Mutect2 on Covered Regions]
    
    Method -->|Traditional| HaplotypeCaller_RNA[HaplotypeCaller<br/>RNA settings]
    
    Mutect2 --> Filter_RNA[Filter RNA Variants]
    HaplotypeCaller_RNA --> Filter_RNA
    
    Filter_RNA --> ExtractPASS[Extract PASS Variants]
    ExtractPASS --> AnnotateRNA[SnpEff Annotation<br/>RNA-aware]
    
    AnnotateRNA --> AddRSID_RNA[Add dbSNP IDs]
    AddRSID_RNA --> AddGenes_RNA[Add Gene Info]
    AddGenes_RNA --> ExportRNA[Export RNA Variants to TSV]
    
    ExportRNA --> EndRNA[End: RNA Variants]
```


# Quality Control & Reporting Flow

```mermaid
graph TD
    StartQC[QC Pipeline] --> RawFastQC[FastQC: Raw FASTQ]
    RawFastQC --> TrimQC[Trimmomatic QC]
    TrimQC --> TrimmedFastQC[FastQC: Trimmed FASTQ]
    
    TrimmedFastQC --> AlignmentQC[Alignment Metrics]
    AlignmentQC --> DuplicateQC[Duplicate Metrics]
    
    DuplicateQC --> VariantQC[Variant Calling Metrics]
    VariantQC --> AnnotationQC[Annotation Statistics]
    
    AnnotationQC --> MultiQC_Report[MultiQC Report Generation]
    MultiQC_Report --> FinalReport[Final QC Report]
    
    FinalReport --> EndQC[End: QC Complete]
```
#  Error Handling & Resume Flow
```mermaid
graph TD
    StartResume[Resume Pipeline] --> ForEachSample[For Each Sample]
    
    ForEachSample --> CheckFiles{Check Output Files}
    
    CheckFiles --> SRA_Exists{SRA exists?}
    SRA_Exists -->|Yes| SkipDownload
    SRA_Exists -->|No| DownloadNeeded
    
    CheckFiles --> FASTQ_Exists{FASTQ exists?}
    FASTQ_Exists -->|Yes| SkipExtract
    FASTQ_Exists -->|No| ExtractNeeded
    
    CheckFiles --> BAM_Exists{BAM exists?}
    BAM_Exists -->|Yes| SkipAlign
    BAM_Exists -->|No| AlignNeeded
    
    CheckFiles --> VCF_Exists{VCF/TSV exists?}
    VCF_Exists -->|Yes| SkipProcess
    VCF_Exists -->|No| ProcessNeeded
    
    SkipDownload --> FASTQ_Exists
    DownloadNeeded --> ExtractNeeded
    SkipExtract --> BAM_Exists
    ExtractNeeded --> AlignNeeded
    SkipAlign --> VCF_Exists
    AlignNeeded --> ProcessNeeded
    SkipProcess --> NextResume{More Samples?}
    ProcessNeeded --> NextResume
    
    NextResume -->|Yes| ForEachSample
    NextResume -->|No| EndResume[End: Resume Complete]
```
# Configuration & Environment Flow
```mermaid
graph TD
    StartConfig[Configuration] --> BaseDirs[Setup Base Directories]
    BaseDirs --> ToolPaths[Set Tool Paths]
    ToolPaths --> RefPaths[Set Reference Paths]
    
    RefPaths --> CheckRefs{References Downloaded?}
    CheckRefs -->|No| DownloadRefs[Download References]
    CheckRefs -->|Yes| SkipRefs[Skip Download]
    
    DownloadRefs --> IndexRefs[Index References]
    SkipRefs --> IndexRefs
    
    IndexRefs --> CreateGeneBed[Create Gene BED File]
    CreateGeneBed --> FuncotatorDS{Funcotator Data Sources?}
    
    FuncotatorDS -->|No| DownloadFunc[Download Funcotator DB]
    FuncotatorDS -->|Yes| SkipFunc[Skip Download]
    
    DownloadFunc --> CreateEnv[Create Environment Script]
    SkipFunc --> CreateEnv
    
    CreateEnv --> EndConfig[End: Configuration Ready]
```
# Differential Analysis- Overview
```mermaid
flowchart TD
    Start([Start MDD2 Differential Analysis]) --> Config[Load Configuration]
    
    Config --> CheckData{Data Check}
    
    CheckData -->|TSV Files Found| LoadMap[Load Project Mapping File]
    CheckData -->|No Files| Error[Error: Run Pipeline First<br>./mdd2.sh run]
    
    LoadMap --> SetParams[Set Analysis Parameters]
    
    subgraph Parameters [Analysis Parameters]
        P1[Threshold: 95%]
        P2[Min Samples: 5]
        P3[Top N: 15]
        P4[Memory Limit: 2GB]
    end
    
    SetParams --> Parameters
    Parameters --> Init[Initialize Processing]
    
    Init --> P1Start[Phase 1 Start]
    Error --> End([Exit])
    
    style Start fill:#e1f5fe
    style Parameters fill:#f3e5f5
    style P1Start fill:#c8e6c9
    style Error fill:#ffcdd2
```
# Differential Analysis - Stream Processing Engine
```mermaid
flowchart TD
    P1Start[Phase 1: Stream Processing] --> GetSamples[Get Sample List]
    
    GetSamples --> ProcessLoop{Process Samples Loop}
    
    ProcessLoop -->|Next Sample| ReadTSV[Read TSV.gz File]
    ReadTSV --> Decompress[Stream Decompress]
    
    Decompress --> Parse[Parse Each Line]
    
    subgraph Parsing [Parse SNP Data]
        P1[Extract CHROM]
        P2[Extract POS]
        P3[Extract REF]
        P4[Extract ALT]
        P5[Extract GENE]
        P6[Create SNP Key<br>CHROM:POS:REF:ALT]
    end
    
    Parse --> Parsing
    
    Parsing --> GroupCheck{Check Group Membership}
    
    GroupCheck -->|Control| ToCtrl[Add to Control Stream]
    GroupCheck -->|MDD| ToMDD[Add to MDD Stream]
    GroupCheck -->|Unknown| Skip[Skip - Log Warning]
    
    ToCtrl --> CountCtrl[Count Control SNPs]
    ToMDD --> CountMDD[Count MDD SNPs]
    Skip --> ProcessLoop
    
    CountCtrl --> StoreCtrl[Store in Control Temp File]
    CountMDD --> StoreMDD[Store in MDD Temp File]
    
    StoreCtrl --> ProcessLoop
    StoreMDD --> ProcessLoop
    
    ProcessLoop -->|All Processed| Merge[Merge Results]
    
    Merge --> Output1[Output: Group-specific Counts]
    
    style P1Start fill:#e1f5fe
    style Parsing fill:#f3e5f5
    style GroupCheck fill:#fff3e0
    style Output1 fill:#c8e6c9
```
# Memory Management
```mermaid
flowchart LR
    Problem[Problem: 1000 samples × 100K SNPs<br>= 100M entries] --> Solution[Solution: Streaming Architecture]
    
    subgraph Streaming [Stream Processing Engine]
        direction TB
        S1[Read Files in Stream] --> S2[Process Line by Line]
        S2 --> S3[Immediate Sorting]
        S3 --> S4[Write to Temp Files]
    end
    
    Solution --> Streaming
    
    subgraph DiskBased [Disk-Based Operations]
        D1[External Sorting<br>sort command] --> D2[Counting<br>uniq -c]
        D2 --> D3[Merging<br>cat + sort]
    end
    
    Streaming --> DiskBased
    
    subgraph GroupWise [Group-Wise Processing]
        G1[Process Control Separately] --> G2[Process MDD Separately]
        G2 --> G3[Merge Only Results]
    end
    
    DiskBased --> GroupWise
    
    subgraph EarlyFilter [Early Filtering]
        F1[Filter by Frequency] --> F2[Keep Only ≥5%]
        F2 --> F3[Discard Rare SNPs]
    end
    
    GroupWise --> EarlyFilter
    
    subgraph TopN [Top-N Storage]
        T1[Store Only Top 15] --> T2[Discard Rest]
        T2 --> T3[Final Size: <1MB]
    end
    
    EarlyFilter --> TopN
    
    TopN --> Result[Result: Efficient Processing]
    
    style Problem fill:#ffebee
    style Result fill:#c8e6c9
    style Streaming fill:#e1f5fe
    style DiskBased fill:#f3e5f5
```
# Differential Analysis Logic
```mermaid
flowchart TD
    StartDiff[Differential Analysis] --> LoadCounts[Load Group Counts]
    
    LoadCounts --> Calculate[Calculate Frequencies]
    
    subgraph FreqCalc [Frequency Calculation]
        FC1[Control Freq = Count / Total Control]
        FC2[MDD Freq = Count / Total MDD]
        FC3[Apply 95% Threshold]
    end
    
    Calculate --> FreqCalc
    
    FreqCalc --> Analyze{Analyze Each SNP}
    
    Analyze -->|Control ≥95% & MDD ≤5%| DiffCtrl[Differential Control SNP]
    Analyze -->|MDD ≥95% & Control ≤5%| DiffMDD[Differential MDD SNP]
    Analyze -->|Other| Common[Common SNP - Discard]
    
    DiffCtrl --> StoreCtrl[Store Control-specific]
    DiffMDD --> StoreMDD[Store MDD-specific]
    
    StoreCtrl --> StatsPrep[Prepare for Statistics]
    StoreMDD --> StatsPrep
    
    StatsPrep --> Output2[Output: Differential SNPs]
    
    style StartDiff fill:#e1f5fe
    style FreqCalc fill:#f3e5f5
    style Analyze fill:#fff3e0
    style DiffCtrl fill:#bbdefb
    style DiffMDD fill:#ffcdd2
    style Output2 fill:#c8e6c9
```
# Statistical Analysis Pipeline
```mermaid
flowchart TD
    StartStats[Statistical Analysis] --> LoadData[Load Differential SNPs]
    
    LoadData --> Annotate[Annotate with Gene Info]
    
    subgraph StatsCalc [Statistical Calculations]
        SC1[2×2 Contingency Table<br>Present/Absent × Control/MDD]
        SC2[Fisher's Exact Test]
        SC3[Calculate P-value]
        SC4[Calculate Odds Ratio]
        SC5[Calculate Confidence Interval]
    end
    
    Annotate --> StatsCalc
    
    StatsCalc --> Rank[Rank by Significance]
    
    Rank --> FilterTop{Filter Top Results}
    
    FilterTop -->|Top 15 Control| TopCtrl[Control Top SNPs]
    FilterTop -->|Top 15 MDD| TopMDD[MDD Top SNPs]
    
    TopCtrl --> FormatCtrl[Format Results]
    TopMDD --> FormatMDD[Format Results]
    
    FormatCtrl --> OutputStats[Statistical Results]
    FormatMDD --> OutputStats
    
    style StartStats fill:#e1f5fe
    style StatsCalc fill:#f3e5f5
    style FilterTop fill:#fff3e0
    style OutputStats fill:#c8e6c9
```
#  Report Generation
```mermaid
flowchart TD
    StartStats[Statistical Analysis] --> LoadData[Load Differential SNPs]
    
    LoadData --> Annotate[Annotate with Gene Info]
    
    subgraph StatsCalc [Statistical Calculations]
        SC1[2×2 Contingency Table<br>Present/Absent × Control/MDD]
        SC2[Fisher's Exact Test]
        SC3[Calculate P-value]
        SC4[Calculate Odds Ratio]
        SC5[Calculate Confidence Interval]
    end
    
    Annotate --> StatsCalc
    
    StatsCalc --> Rank[Rank by Significance]
    
    Rank --> FilterTop{Filter Top Results}
    
    FilterTop -->|Top 15 Control| TopCtrl[Control Top SNPs]
    FilterTop -->|Top 15 MDD| TopMDD[MDD Top SNPs]
    
    TopCtrl --> FormatCtrl[Format Results]
    TopMDD --> FormatMDD[Format Results]
    
    FormatCtrl --> OutputStats[Statistical Results]
    FormatMDD --> OutputStats
    
    style StartStats fill:#e1f5fe
    style StatsCalc fill:#f3e5f5
    style FilterTop fill:#fff3e0
    style OutputStats fill:#c8e6c9
```
# Error Handling & Recovery
```mermaid
flowchart TD
    StartError[Error Handling System] --> TryProcess[Attempt Processing]
    
    TryProcess --> ErrorCheck{Error Occurs?}
    
    ErrorCheck -->|No Error| Continue[Continue Normally]
    
    ErrorCheck -->|File Not Found| HandleFile[Handle Missing File]
    ErrorCheck -->|Memory Error| HandleMem[Handle Memory Issue]
    ErrorCheck -->|Parse Error| HandleParse[Handle Parse Error]
    ErrorCheck -->|Unknown Error| HandleUnknown[Handle Unknown Error]
    
    HandleFile --> Log1[Log Warning] --> Skip1[Skip Sample] --> Recover1[Continue with Rest]
    HandleMem --> Log2[Log Error] --> Adjust1[Reduce Batch Size] --> Retry1[Retry]
    HandleParse --> Log3[Log Error] --> Skip2[Skip Line] --> Recover2[Continue File]
    HandleUnknown --> Log4[Log Critical] --> SavePartial[Save Partial Results] --> ExitGraceful[Exit Gracefully]
    
    Recover1 --> CheckProgress{Enough Data?}
    Retry1 --> CheckProgress
    Recover2 --> CheckProgress
    ExitGraceful --> EndError([Exit with Error])
    
    CheckProgress -->|Yes ≥Min Samples| ContinueAnalysis[Continue Analysis]
    CheckProgress -->|No| Abort[Abort Analysis]
    
    Continue --> ContinueAnalysis
    ContinueAnalysis --> GeneratePartial[Generate Partial Report]
    
    GeneratePartial --> EndPartial([Analysis Complete - Partial Results])
    Abort --> EndError
    
    style StartError fill:#e1f5fe
    style ErrorCheck fill:#fff3e0
    style HandleFile fill:#ffebee
    style HandleMem fill:#ffebee
    style HandleParse fill:#ffebee
    style HandleUnknown fill:#ffebee
    style ContinueAnalysis fill:#c8e6c9
    style GeneratePartial fill:#c8e6c9
```
# Data flow diagram
```mermaid
flowchart TD
    Input[TSV.gz Files<br>100s of samples] --> Decompress[Stream Decompress]
    
    Decompress --> Parse[Parse Each Line]
    Parse --> Extract[Extract: CHROM, POS, REF, ALT, GENE]
    
    Extract --> CreateKey[Create SNP Key<br>CHROM:POS:REF:ALT]
    
    subgraph GroupSeparation [Group Separation]
        GS1{Control or MDD?}
        GS1 -->|Control| C1[Add to Control Stream]
        GS1 -->|MDD| M1[Add to MDD Stream]
    end
    
    CreateKey --> GroupSeparation
    
    C1 --> CountCtrl[Count per SNP in Control]
    M1 --> CountMDD[Count per SNP in MDD]
    
    CountCtrl --> MergeCtrl[Merge All Control Counts]
    CountMDD --> MergeMDD[Merge All MDD Counts]
    
    MergeCtrl --> FreqCtrl[Calculate Frequency<br>Control]
    MergeMDD --> FreqMDD[Calculate Frequency<br>MDD]
    
    FreqCtrl --> Compare{Compare Frequencies}
    FreqMDD --> Compare
    
    Compare -->|Ctrl ≥95%, MDD ≤5%| DiffCtrl[Differential: Control]
    Compare -->|MDD ≥95%, Ctrl ≤5%| DiffMDD[Differential: MDD]
    Compare -->|Other| Common[Common - Discard]
    
    DiffCtrl --> StatsCtrl[Statistical Analysis<br>Control-specific]
    DiffMDD --> StatsMDD[Statistical Analysis<br>MDD-specific]
    
    StatsCtrl --> RankCtrl[Rank by Significance]
    StatsMDD --> RankMDD[Rank by Significance]
    
    RankCtrl --> TopCtrl[Top 15 Control SNPs]
    RankMDD --> TopMDD[Top 15 MDD SNPs]
    
    TopCtrl --> ReportCtrl[Report Generation]
    TopMDD --> ReportMDD[Report Generation]
    
    ReportCtrl --> FinalExcel[Final Excel Report]
    ReportMDD --> FinalExcel
    
    FinalExcel --> Output[Output: Small, Focused Results<br><100KB vs GB matrix]
    
    style Input fill:#e1f5fe
    style Output fill:#c8e6c9
    style DiffCtrl fill:#bbdefb
    style DiffMDD fill:#ffcdd2
```
# Error recovery code
```mermaid
flowchart TD
    ProcessSample[Process Sample] --> CheckFile{File Exists?}
    CheckFile -->|Yes| ProcessOK[Process Normally]
    CheckFile -->|No| LogError[Log Error]
    
    ProcessOK --> Continue[Continue to Next]
    LogError --> Skip[Skip This Sample]
    
    Skip --> CheckProgress{Enough Samples Processed?}
    CheckProgress -->|Yes, ≥Min Required| ContinueAnalysis[Continue Analysis]
    CheckProgress -->|No| Abort[Abort Analysis]
    
    Continue --> ContinueAnalysis
    ContinueAnalysis --> Final[Generate Partial Results]
    Abort --> Cleanup[Cleanup & Exit]
    
    style ProcessOK fill:#c8e6c9
    style LogError fill:#ffebee
    style Final fill:#fff3e0
```
