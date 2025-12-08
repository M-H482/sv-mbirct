#ifndef _ACOMP_STREAM_H_
#define _ACOMP_STREAM_H_

#include "MBIRModularDefs.h"
#include "A_comp.h" // 复用原有的结构体定义

/* * 流式计算主函数
 * * memory_limit_gb: 允许使用的最大内存 (GB)。例如 1024.0 或 64.0
 * verboseLevel: 输出详细程度
 */
void AmatrixComputeToFile_Stream(
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    char *Amatrix_fname,
    char verboseLevel,
    double memory_limit_gb);

#endif
