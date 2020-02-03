// (c) 2008-2009 Miaoxin Li
// This file is distributed as part of the IGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.

// Permission is granted for you to use this file to compile IGG.

// All computer programs have bugs. Use this file at your own risk.
// Saturday, January 17, 2009
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.eqtlsimulator;

/**
 *
 * @author mxli
 */
public interface Constants {

    public static final String JAVA_VERSION = System.getProperty("java.version");
    public static final String IGG_VERSION = "3.0 ; Java/" + JAVA_VERSION;
    public static final String AUTHOR_STRING = "Miao-Xin Li; limx54@yahoo.com";
    public static final String IGG_WEBSITE_STRING = "http://bioinfo.hku.hk/iggweb";
    public static final String IGG_CITATION_STRING = "<html>Li et al. (2007) IGG: A Tool to Integrate Gene chips for Genetic Studies<br>" +
            " Bioinformatics 23(22):3105-3107</html>";
    public static final String RESOURCE_PATH = "resources/";
    public final String IGG_RESOURCE_CHIP_FOLDER = "resources/chip_annotation/";
    public final String IGG_RESOURCE_HAPMAP_FOLDER = "resources/hapmap/";
    public static final String IGG_RESOURCE_URL = "http://bioinfo.hku.hk:13080/iggweb/download/";
    static int FILE_READER_BUFFER_SIZE = 10 * 1024 * 1024;
    static int HAP_GENOTYPE_LOAD_FACTOR = 10 * 1024;
    static int CHIP_GENOTYPE_LOAD_FACTOR = 1024;
    static final String HELP_OUTPUT = " Command line options\n" +
            "-h, -help                       Print this message\n" +
            "-memory <memsize>               allocates <memsize> megabytes of memory (default 512)\n" +
            "-n, -nogui                      Command line output only\n" +
            "-q, -quiet                      Quiet mode- doesnt print any warnings or information to screen\n" +
            "-log <filename>                 Specify a logfile name (defaults to haploview.log if no name specified)\n";
    public static final String[] ILLUMINA_CHIPS = {"Illumina_HumanHap1M", "Illumina_HumanHap650Y", "Illumina_HumanHap610-Quad",
        "Illumina_HumanHap550-Duo", "Illumina_HumanHapCNV370", "Illumina_HumanHap300-Duo",
        "Illumina_Human_Linkage_IVb_Panel", "Illumina_HumanLinkage-12"
    };
    public static final String[] AFFY_CHIPS = {"Affymetrix_GenomeWide_Human_SNP_Array_6.0", "Affymetrix_GenomeWide_Human_SNP_Array_5.0",
        "Affymetrix_Mapping_250K_Nsp_Array", "Affymetrix_Mapping_250K_Sty_Array", "Affymetrix_Mapping_50K_Hind_240_Array",
        "Affymetrix_Mapping_50K_Xba_240_Array", "Affymetrix_Mapping_10K_2.0_Array"
    };
    public static final String[] ALL_CHIP_InSizeOrder = {"Illumina_HumanHap1M", "Affymetrix_GenomeWide_Human_SNP_Array_6.0", "Illumina_HumanHap650Y", "Illumina_HumanHap610-Quad",
        "Illumina_HumanHap550-Duo", "Affymetrix_GenomeWide_Human_SNP_Array_5.0",
        "Illumina_HumanHapCNV370", "Illumina_HumanHap300-Duo", "Affymetrix_Mapping_250K_Nsp_Array", "Affymetrix_Mapping_250K_Sty_Array", "Affymetrix_Mapping_50K_Hind_240_Array",
        "Affymetrix_Mapping_50K_Xba_240_Array", "Affymetrix_Mapping_10K_2.0_Array", "Illumina_Human_Linkage_IVb_Panel", "Illumina_HumanLinkage-12"
    };
    public static final String[] CHROM_NAME = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
        "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT", "XY"
    };
    public static final int[] STANDARD_CHROM_NUM = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    public static final String[] HAPMAP_POPULATIONS = {"_CHB", "_JPT", "_YRI", "_CEU", "_ASW", "_CHD", "_GIH", "_LWK", "_MEX", "_MKK", "_TSI", "JPT+CHB"};
    public static final String[] HAPMAP_ANNO = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
        "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"
    };
    public static int DEFAULT_BUF_SIZE = 2500;
    public static int DEFAULT_GTY_BLOCKNUM = 4;
    public static char MISSING_ALLELE_NAME = 'X';
    public static char MISSING_STRAND_NAME = '0';
    public static char DEFAULT_MISSING_GTY_NAME = '0';
    public static double DEFAULT_MISSING_DOUBLE=-9.99;
    //download setting
    public static int MAX_TOTAL_CONNECTIONS = 1000;
    public static int MAX_PER_ROUTE_CONNECTIONS = 50;
    public static int MAX_THREAD_NUM = 8;
    public static int FILE_SEGEMENT_NUM = 10;

}
