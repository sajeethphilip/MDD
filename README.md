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
    A[Start MDD2 Analysis] --> B[Discover Project Structure]
    B --> C{Multi-directory Project?}
    C -->|Yes| D[Setup Multi-directory Analysis]
    C -->|No| E[Setup Single Directory Analysis]
    
    D --> F[Select MDD Directory]
    D --> G[Select Control Directory]
    F & G --> H[Collect & Unify TSV Files]
    
    E --> I[Collect TSV Files from Single Directory]
    
    H --> J[Create Master Directory Structure]
    I --> J
    
    J --> K[Create Project Mapping]
    K --> L[Run Differential SNP Analysis]
    
    L --> M{Success?}
    M -->|Yes| N[Generate Comprehensive Reports]
    M -->|No| O[Error Handling & Logging]
    
    N --> P[Output Results]
    O --> Q[Display Error & Exit]
    
    P --> R[End: Analysis Complete]
    Q --> S[End: Analysis Failed]
  
```
