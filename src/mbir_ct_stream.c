#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/time.h>

#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "initialize.h"
#include "A_comp_stream.h"

/* 简化的命令行参数结构体，仅包含矩阵生成所需字段 */
struct CmdLine {
    char ImageParamsFile[1024];
    char SinoParamsFile[1024];
    char SysMatrixFile[1024];
    int ImageParamsFileFlag;
    int SinoParamsFileFlag;
    int SysMatrixFileFlag;
    int verboseLevel;
    double memory_limit_gb;
};

void readCmdLine(int argc, char *argv[], struct CmdLine *cmdline);
void printCmdLineUsage(char *ExecFileName);

int main(int argc, char *argv[])
{
    struct CmdLine cmdline;
    struct ImageParams3D imgparams;
    struct SinoParams3DParallel sinoparams;
    char fname[1064];
    struct timeval tm0, tm2;
    unsigned long long tdiff;

    gettimeofday(&tm0, NULL);

    /* 1. 解析命令行参数 */
    readCmdLine(argc, argv, &cmdline);

    if(cmdline.verboseLevel) {
        fprintf(stdout, "SUPER-VOXEL MBIR CT - SYSTEM MATRIX GENERATOR (STREAMING)\n");
        fprintf(stdout, "Build time: %s, %s\n", __DATE__,  __TIME__);
    }

    /* 2. 读取参数文件 */
    if(cmdline.verboseLevel) fprintf(stdout, "Reading parameter files...\n");
    ReadSinoParams3DParallel(cmdline.SinoParamsFile, &sinoparams);
    ReadImageParams3D(cmdline.ImageParamsFile, &imgparams);

    /* 3. 执行流式矩阵计算 */
    if(cmdline.SysMatrixFileFlag)
    {
        sprintf(fname, "%s.2Dsvmatrix", cmdline.SysMatrixFile);
        
        fprintf(stdout, "------------------------------------------------\n");
        fprintf(stdout, "Operation: Compute System Matrix (Streaming)\n");
        fprintf(stdout, "Output File: %s\n", fname);
        fprintf(stdout, "Memory Limit: %.2f GB\n", cmdline.memory_limit_gb);
        fprintf(stdout, "------------------------------------------------\n");

        AmatrixComputeToFile_Stream(
            imgparams, 
            sinoparams, 
            fname, 
            cmdline.verboseLevel, 
            cmdline.memory_limit_gb
        );
    }
    else {
        /* 理论上 readCmdLine 已经检查了，这里是双重保险 */
        fprintf(stderr, "Error: System matrix output file not specified (-m).\n");
        exit(-1);
    }

    if(cmdline.verboseLevel) {
        gettimeofday(&tm2, NULL);
        tdiff = 1000 * (tm2.tv_sec - tm0.tv_sec) + (tm2.tv_usec - tm0.tv_usec) / 1000;
        fprintf(stdout, "Done. Total run time = %llu ms\n", tdiff);
    }

    return 0;
}

void readCmdLine(int argc, char *argv[], struct CmdLine *cmdline)
{
    char ch;
    
    /* 设置默认值 */
    cmdline->SinoParamsFileFlag = 0;
    cmdline->ImageParamsFileFlag = 0;
    cmdline->SysMatrixFileFlag = 0;
    cmdline->verboseLevel = 1;
    cmdline->memory_limit_gb = 1024; /* 默认 1024GB */

    if(argc == 1) {
        printCmdLineUsage(argv[0]);
        exit(0);
    }
    
    /* 解析参数: 
       i: imgparams
       j: sinoparams
       m: output matrix base name
       v: verbose
       M: Memory limit (GB) - 新增
    */
    while ((ch = getopt(argc, argv, "i:j:m:v:M:h")) != EOF)
    {
        switch (ch)
        {
            case 'i':
                cmdline->ImageParamsFileFlag = 1;
                sprintf(cmdline->ImageParamsFile, "%s", optarg);
                break;
            case 'j':
                cmdline->SinoParamsFileFlag = 1;
                sprintf(cmdline->SinoParamsFile, "%s", optarg);
                break;
            case 'm':
                cmdline->SysMatrixFileFlag = 1;
                sprintf(cmdline->SysMatrixFile, "%s", optarg);
                break;
            case 'v':
                sscanf(optarg, "%d", &cmdline->verboseLevel);
                break;
            case 'M':
                sscanf(optarg, "%lf", &cmdline->memory_limit_gb);
                break;
            case 'h':
                printCmdLineUsage(argv[0]);
                exit(0);
                break;
            default:
                fprintf(stderr, "Unknown option: -%c\n", ch);
                printCmdLineUsage(argv[0]);
                exit(-1);
        }
    }

    /* 检查必要参数 */
    if(!cmdline->SinoParamsFileFlag || !cmdline->ImageParamsFileFlag || !cmdline->SysMatrixFileFlag){
        fprintf(stderr, "Error: Missing mandatory arguments.\n");
        fprintf(stderr, "You must specify: -i <imgparams> -j <sinoparams> -m <output_matrix>\n");
        exit(-1);
    }
}

void printCmdLineUsage(char *ExecFileName)
{
    fprintf(stdout, "\nUsage: %s -i <imgparams> -j <sinoparams> -m <output_matrix_base> [options]\n", ExecFileName);
    fprintf(stdout, "Options:\n");
    fprintf(stdout, "  -M <float>  Set memory limit in GB (Default: 64.0)\n");
    fprintf(stdout, "  -v <int>    Set verbosity level (0: quiet, 1: default, 2: verbose)\n");
    fprintf(stdout, "Example:\n");
    fprintf(stdout, "  %s -i data.imgparams -j data.sinoparams -m matrix_out -M 128.0\n\n", ExecFileName);
}
