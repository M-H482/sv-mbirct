#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "MBIRModularDefs.h"
#include "MBIRModularUtils.h"
#include "allocate.h"
#include "initialize.h"
#include "A_comp.h"
#include "A_comp_stream.h"

#ifndef MSVC
#include <sys/time.h>
#endif

#define PI 3.1415926535897932384
#define LEN_DET 101
#define LEN_PIX 511
#define LEN_ANG 511

/* --- 内部 Helper 函数 (保持不变) --- */

double angle_mod_s(double theta, double lower, double upper)
{
    double x,z,interval;
    interval = upper-lower;
    x = theta-lower;
    z = fmod(x,interval);
    if(x<0)
        z += interval;
    return z + lower;
}

float UnitPixelProj_s(float angle, float t)
{
    float radius = 0.7071067811865476;
    float proj,d1,d2,tmag;
    angle = angle_mod_s(angle,-PI/4.0,PI/4.0);
    if(angle>0) {
        d1 = radius*cosf(angle + PI/4.0);
        d2 = radius*cosf(angle - PI/4.0);
    } else {
        d1 = radius*cosf(angle - PI/4.0);
        d2 = radius*cosf(angle + PI/4.0);
    }
    tmag = fabs(t);
    if(tmag >= d2) proj = 0.0;
    else if (tmag <= d1) proj = 1.0/cosf(angle);
    else proj = 1.0/cosf(angle) * (d2-tmag)/(d2-d1);
    return proj;
}

float **ComputePixelProfLookup_s(float Deltaxy)
{
    int i,j;
    float angle,t;
    float **pix_prof;
    pix_prof = (float **)get_img(LEN_PIX, LEN_ANG, sizeof(float));
    for (i = 0; i < LEN_ANG; i++) {
        angle = (PI/2.0)/LEN_ANG * i;
        for (j = 0; j < LEN_PIX; j++) {
            t = 1.0/LEN_PIX * j;
            pix_prof[i][j] = Deltaxy * UnitPixelProj_s(angle, t);
        }
    }
    return pix_prof;
}

float PixProjLookup_s(float **pix_prof, float Deltaxy, float angle, float t)
{
    int i_ang, i_t;
    float proj;
    angle = angle_mod_s(angle,0.0,PI/2.0);
    i_ang = (int) (angle/(PI/2.0)*LEN_ANG + 0.5);
    if(i_ang == LEN_ANG) i_ang = 0;
    i_t = (int) (fabs(t)/Deltaxy*LEN_PIX + 0.5);
    if(i_t >= LEN_PIX) proj = 0.0;
    else proj = pix_prof[i_ang][i_t];
    return proj;
}

void A_comp_ij_s(
    int im_row,
    int im_col,
    struct SinoParams3DParallel *sinoparams,
    struct ImageParams3D *imgparams,
    float **pix_prof,
    struct ACol *A_col,float *A_Values)
{
    static int first_call=1;
    static float dprof[LEN_DET];
    int i, k, pr, ind_min, ind_max, proj_count;
    float t_0, x_0, y_0, x, y;
    float t, t_pix=0.0, t_min, t_max, t_start;
    float Aval, detSampleD;
    float r_sd, r_si=1.0, x_s, y_s, theta=0.0, alpha=0.0, D=1.0, M=1.0;

    float Deltaxy = imgparams->Deltaxy;
    int NChannels = sinoparams->NChannels;
    float DeltaChannel;

    if (first_call == 1) {
        first_call = 0;
        for (k = 0; k < LEN_DET; k++) dprof[k] = 1.0/(LEN_DET);
    }

    r_sd = sinoparams->DistSourceDetector;
    r_si = r_sd / sinoparams->Magnification;
    DeltaChannel = sinoparams->DeltaChannel;
    if(sinoparams->Geometry == 1) DeltaChannel = sinoparams->DeltaChannel / r_sd;

    if(LEN_DET > 1) detSampleD = DeltaChannel/(LEN_DET-1);
    else detSampleD = DeltaChannel/2.0;

    t_0 = -(NChannels-1)*DeltaChannel/2.0 - sinoparams->CenterOffset * DeltaChannel;
    x_0 = -(imgparams->Nx-1)*Deltaxy/2.0;
    y_0 = -(imgparams->Ny-1)*Deltaxy/2.0;

    y = y_0 + im_row*Deltaxy;
    x = x_0 + im_col*Deltaxy;

    proj_count = 0;
    for (pr = 0; pr < sinoparams->NViews; pr++)
    {
        int countTemp=proj_count;
        int write=1;
        int minCount=0;
        float view_angle = sinoparams->ViewAngles[pr];

        if(sinoparams->Geometry == 1) { // curved
            x_s = r_si * cosf(view_angle);
            y_s = r_si * sinf(view_angle);
            theta = atan2(y_s-y, x_s-x);
            alpha = angle_mod_s(theta - view_angle,-PI,PI);
            D = sqrt((x_s-x)*(x_s-x) + (y_s-y)*(y_s-y));
            t_pix = alpha;
            t_min = t_pix - Deltaxy/D;
            t_max = t_pix + Deltaxy/D;
        } else if(sinoparams->Geometry == 2) { // flat
            x_s = r_si * cosf(view_angle);
            y_s = r_si * sinf(view_angle);
            theta = atan2(y_s-y, x_s-x);
            alpha = angle_mod_s(theta - view_angle,-PI,PI);
            D = sqrt((x_s-x)*(x_s-x) + (y_s-y)*(y_s-y));
            M = r_sd/cosf(alpha) / D;
            t_pix = r_sd*tanf(alpha);
            t_min = t_pix - Deltaxy*M;
            t_max = t_pix + Deltaxy*M;
        } else { // parallel
            t_pix = y*cosf(view_angle) - x*sinf(view_angle);
            t_min = t_pix - Deltaxy;
            t_max = t_pix + Deltaxy;
        }

        ind_min = ceil((t_min-t_0)/DeltaChannel - 0.5);
        ind_max= floor((t_max-t_0)/DeltaChannel + 0.5);

        if(ind_max<0 || ind_min>NChannels-1) {
            A_col->countTheta[pr]=0;
            A_col->minIndex[pr]=0;
            continue;
        }

        ind_min = (ind_min<0) ? 0 : ind_min;
        ind_max = (ind_max>=NChannels) ? NChannels-1 : ind_max;

        for (i = ind_min; i <= ind_max; i++) {
            Aval = 0;
            t_start = t_0 - DeltaChannel/2.0 + i*DeltaChannel;
            for (k = 0; k < LEN_DET; k++) {
                t = t_start + k*detSampleD;
                if(sinoparams->Geometry == 1)
                    Aval += dprof[k]*PixProjLookup_s(pix_prof, Deltaxy, theta, (t-alpha)*D);
                else if(sinoparams->Geometry == 2)
                    Aval += dprof[k]*PixProjLookup_s(pix_prof, Deltaxy, theta, (t-t_pix)*cosf(alpha)/M);
                else
                    Aval += dprof[k]*PixProjLookup_s(pix_prof, Deltaxy, view_angle, t-t_pix);
            }
            if (Aval > 0.0) {
                if(write==1) { minCount=i; write=0; }
                A_Values[proj_count] = Aval;
                proj_count++;
            }
        }
        A_col->countTheta[pr] = proj_count-countTemp;
        A_col->minIndex[pr] = minCount;
    }
    A_col->n_index = proj_count;
}

int estimate_batch_size(double mem_gb, struct SVParams svpar, struct SinoParams3DParallel *sinoparams) {
    long long bytes_available = (long long)(mem_gb * 1024 * 1024 * 1024);
    
    int sv_width = 2 * svpar.SVLength + 1;
    int pixels_per_sv = sv_width * sv_width;
    
    long long mem_per_pixel_acol = sinoparams->NViews * (sizeof(chanwidth_t) + sizeof(channel_t));
    long long mem_per_pixel_aval = (long long)(sinoparams->NViews * sinoparams->NChannels * 0.1); 
    if (mem_per_pixel_aval < 1024) mem_per_pixel_aval = 1024;

    long long mem_per_sv = pixels_per_sv * (mem_per_pixel_acol + mem_per_pixel_aval + sizeof(struct ACol) + sizeof(struct AValues_char));
    
    long long safe_mem = (long long)(bytes_available * 0.7);
    
    int batch_size = (int)(safe_mem / mem_per_sv);
    if (batch_size < 1) batch_size = 1;
    if (batch_size > svpar.Nsv) batch_size = svpar.Nsv;
    
    return batch_size;
}

void AmatrixComputeToFile_Stream(
    struct ImageParams3D imgparams,
    struct SinoParams3DParallel sinoparams,
    char *Amatrix_fname,
    char verboseLevel,
    double memory_limit_gb)
{
    struct SVParams svpar;
    FILE *fp;
    int i, j, t;
    
    if(verboseLevel) fprintf(stdout, "Initializing SV Params...\n");
    initSVParams(&svpar, imgparams, sinoparams);

    int Nx = imgparams.Nx;
    int Ny = imgparams.Ny;
    int Nsv = svpar.Nsv;
    int SVLength = svpar.SVLength;
    int NViews = sinoparams.NViews;
    int NViewSets = NViews / svpar.pieceLength;

    for(i=0; i<Nsv; i++) {
        svpar.bandMinMap[i].bandMin = (channel_t*)get_spc(NViews, sizeof(channel_t));
        svpar.bandMaxMap[i].bandMax = (channel_t*)get_spc(NViews, sizeof(channel_t));
        for(t=0; t<NViews; t++) {
            svpar.bandMinMap[i].bandMin[t] = sinoparams.NChannels;
            svpar.bandMaxMap[i].bandMax[t] = 0;
        }
    }

    float *Aval_max_ptr = (float *) get_spc(Nx*Ny, sizeof(float));
    memset(Aval_max_ptr, 0, Nx*Ny*sizeof(float));

    int *order = (int *) mget_spc(Nsv, sizeof(int));
    t=0;
    for(i=0; i<Ny; i+=(SVLength*2-svpar.overlap))
        for(j=0; j<Nx; j+=(SVLength*2-svpar.overlap)) {
            order[t] = i*Nx+j; 
            t++;
        }

    if ((fp = fopen(Amatrix_fname, "wb")) == NULL) {
        fprintf(stderr, "ERROR: can't open file %s for writing.\n", Amatrix_fname);
        exit(-1);
    }

    float **pix_prof = ComputePixelProfLookup_s(imgparams.Deltaxy);
    int batch_size = estimate_batch_size(memory_limit_gb, svpar, &sinoparams);
    
    if(verboseLevel) {
        fprintf(stdout, "Total SVs: %d\n", Nsv);
        fprintf(stdout, "Memory Limit: %.2f GB\n", memory_limit_gb);
        fprintf(stdout, "Batch Size: %d SVs\n", batch_size);
    }

    for (int batch_start = 0; batch_start < Nsv; batch_start += batch_size) 
    {
        int batch_end = batch_start + batch_size;
        if (batch_end > Nsv) batch_end = Nsv;
        int current_batch_count = batch_end - batch_start;

        if(verboseLevel) fprintf(stdout, "Processing Batch SVs %d to %d...\n", batch_start, batch_end-1);

        int y_min = Ny, y_max = 0;
        int *batch_sv_indices = (int*)malloc(current_batch_count * sizeof(int));

        for(int k=0; k<current_batch_count; k++) {
            int sv_idx = batch_start + k;
            int pixel_idx = order[sv_idx];
            int jy = pixel_idx / Nx;
            
            batch_sv_indices[k] = sv_idx;
            if (jy < y_min) y_min = jy;
            if (jy + 2*SVLength + 1 > y_max) y_max = jy + 2*SVLength + 1;
        }
        if (y_max > Ny) y_max = Ny;

        struct ACol **ACol_arr = (struct ACol **)calloc(Ny, sizeof(struct ACol*));
        struct AValues_char **AVal_arr = (struct AValues_char **)calloc(Ny, sizeof(struct AValues_char*));

        for(i = y_min; i < y_max; i++) {
            ACol_arr[i] = (struct ACol *)calloc(Nx, sizeof(struct ACol));
            AVal_arr[i] = (struct AValues_char *)calloc(Nx, sizeof(struct AValues_char));
        }

        #pragma omp parallel private(j)
        {
            struct ACol A_col_sgl;
            A_col_sgl.countTheta = (chanwidth_t *)get_spc(NViews,sizeof(chanwidth_t));
            A_col_sgl.minIndex = (channel_t *)get_spc(NViews,sizeof(channel_t));
            float *A_val_sgl = (float *)get_spc(NViews*sinoparams.NChannels, sizeof(float));
            int r;

            #pragma omp for schedule(dynamic)
            for (i = y_min; i < y_max; i++) {
                for (j = 0; j < Nx; j++) {
                    A_comp_ij_s(i, j, &sinoparams, &imgparams, pix_prof, &A_col_sgl, A_val_sgl);
                    
                    ACol_arr[i][j].n_index = A_col_sgl.n_index;
                    ACol_arr[i][j].countTheta = (chanwidth_t *) get_spc(NViews, sizeof(chanwidth_t));
                    ACol_arr[i][j].minIndex = (channel_t *) get_spc(NViews, sizeof(channel_t));
                    AVal_arr[i][j].val = (unsigned char *) get_spc(A_col_sgl.n_index, sizeof(unsigned char));

                    float maxval = 0.0;
                    if (A_col_sgl.n_index > 0) maxval = A_val_sgl[0];
                    for (r = 0; r < A_col_sgl.n_index; r++) {
                        if(A_val_sgl[r] > maxval) maxval = A_val_sgl[r];
                    }
                    
                    Aval_max_ptr[i*Nx+j] = maxval;

                    /* Fix 1: 量化计算必须严格匹配原始 A_comp.c。
                       原代码: val = (unsigned char)((float/float)*255 + 0.5) 
                       注意原始代码用的是整数 255，会发生类型提升。
                       为了保证二进制完全一致，这里改回 255 (int)，而不是 255.0 (double)。
                    */
                    if (maxval > 0) {
                        for (r=0; r < A_col_sgl.n_index; r++)
                            AVal_arr[i][j].val[r] = (unsigned char)((A_val_sgl[r])/maxval*255 + 0.5);
                    }

                    for (r=0; r < NViews; r++) {
                        ACol_arr[i][j].countTheta[r] = A_col_sgl.countTheta[r];
                        ACol_arr[i][j].minIndex[r] = A_col_sgl.minIndex[r];
                    }
                }
            }
            free((void *)A_val_sgl);
            free((void *)A_col_sgl.countTheta);
            free((void *)A_col_sgl.minIndex);
        }

        /* Fix 2: 关键修复！
           原始 A_piecewise 代码中有一个"Gap Filling"步骤，会填补 countTheta==0 时的 minIndex。
           流式版本必须复刻这个逻辑，否则 bandMin 计算会偏小，导致 Header 不一致。
        */
        #pragma omp parallel for private(j)
        for(i = y_min; i < y_max; i++) {
            for(j = 0; j < Nx; j++) {
                if(ACol_arr[i][j].n_index > 0) {
                    for(int p=0; p<NViews; p++) {
                        if(ACol_arr[i][j].minIndex[p]==0 && ACol_arr[i][j].countTheta[p]==0) {
                            if(p!=0) {
                                ACol_arr[i][j].minIndex[p] = ACol_arr[i][j].minIndex[p-1];
                            } else {
                                int temp_t=0;
                                while(ACol_arr[i][j].minIndex[temp_t] == 0 && temp_t < NViews-1)
                                    temp_t++;
                                ACol_arr[i][j].minIndex[p] = ACol_arr[i][j].minIndex[temp_t];
                            }
                        }
                    }
                }
            }
        }

        struct AValues_char **Batch_SV_Data = (struct AValues_char **)malloc(current_batch_count * sizeof(struct AValues_char*));
        
        #pragma omp parallel for schedule(dynamic)
        for(int k=0; k<current_batch_count; k++) 
        {
            int sv_idx = batch_sv_indices[k];
            int pixel_start_idx = order[sv_idx];
            int jy = pixel_start_idx / Nx;
            int jx = pixel_start_idx % Nx;
            
            int SV_Vol_Count = (2*SVLength+1)*(2*SVLength+1);
            Batch_SV_Data[k] = (struct AValues_char *)calloc(SV_Vol_Count, sizeof(struct AValues_char));

            channel_t *bandMin = svpar.bandMinMap[sv_idx].bandMin;
            
            int SVSize = 0;
            int jx_list[SV_Vol_Count];
            int jy_list[SV_Vol_Count];
            
            for(int jy_new=jy; jy_new<=(jy+2*SVLength); jy_new++)
            for(int jx_new=jx; jx_new<=(jx+2*SVLength); jx_new++) {
                if(jy_new<Ny && jx_new<Nx && ACol_arr[jy_new][jx_new].n_index > 0) {
                     jy_list[SVSize] = jy_new;
                     jx_list[SVSize] = jx_new;
                     SVSize++;
                }
            }

            for(int p=0; p<NViews; p++) bandMin[p] = sinoparams.NChannels;
            channel_t localBandMax[NViews];
            for(int p=0; p<NViews; p++) localBandMax[p] = 0;

            for(int i=0; i<SVSize; i++) {
                int r_y = jy_list[i];
                int r_x = jx_list[i];
                for(int p=0; p<NViews; p++) {
                    if(ACol_arr[r_y][r_x].minIndex[p] < bandMin[p]) bandMin[p] = ACol_arr[r_y][r_x].minIndex[p];
                }
            }
            
            /* Fix 3: 确保 BandMax 初始化逻辑与 BandMin 对齐 */
            for(int p=0; p<NViews; p++) localBandMax[p] = bandMin[p];
            
            for(int i=0; i<SVSize; i++) {
                int r_y = jy_list[i];
                int r_x = jx_list[i];
                for(int p=0; p<NViews; p++) {
                    int end_v = ACol_arr[r_y][r_x].minIndex[p] + ACol_arr[r_y][r_x].countTheta[p];
                    if(end_v > localBandMax[p]) localBandMax[p] = end_v;
                }
            }
            
            memcpy(svpar.bandMaxMap[sv_idx].bandMax, localBandMax, NViews*sizeof(channel_t));

            channel_t bandWidthPW[NViewSets];
            for(int p=0; p<NViewSets; p++) {
                int max_w = 0;
                for(int t=0; t<svpar.pieceLength; t++) {
                     int idx = p*svpar.pieceLength+t;
                     if(idx < NViews) {
                        int w = localBandMax[idx] - bandMin[idx];
                        if(w > max_w) max_w = w;
                     }
                }
                bandWidthPW[p] = max_w;
            }
            
            for(int p=0; p<NViews; p++) {
                if((bandMin[p] + bandWidthPW[p/svpar.pieceLength]) >= sinoparams.NChannels)
                    bandMin[p] = sinoparams.NChannels - bandWidthPW[p/svpar.pieceLength];
            }

            for(int dy=0; dy <= 2*SVLength; dy++) {
                for(int dx=0; dx <= 2*SVLength; dx++) {
                    int v_idx = dy*(2*SVLength+1) + dx;
                    int r_y = jy + dy;
                    int r_x = jx + dx;
                    
                    if(r_y >= Ny || r_x >= Nx || ACol_arr[r_y][r_x].n_index == 0) {
                        Batch_SV_Data[k][v_idx].length = 0;
                        continue;
                    }

                    int totalLen = 0;
                    channel_t *pwMin = (channel_t*)malloc(NViewSets * sizeof(channel_t));
                    channel_t *pwWidth = (channel_t*)malloc(NViewSets * sizeof(channel_t));
                    
                    for(int p=0; p<NViewSets; p++) {
                        int p_start = p * svpar.pieceLength;
                        int min_v = ACol_arr[r_y][r_x].minIndex[p_start] - bandMin[p_start];
                        int max_v = min_v + ACol_arr[r_y][r_x].countTheta[p_start];
                        
                        for(int t=0; t<svpar.pieceLength; t++) {
                            int idx = p_start + t;
                            if(idx < NViews) {
                                int v0 = ACol_arr[r_y][r_x].minIndex[idx] - bandMin[idx];
                                int v1 = v0 + ACol_arr[r_y][r_x].countTheta[idx];
                                if(v0 < min_v) min_v = v0;
                                if(v1 > max_v) max_v = v1;
                            }
                        }
                        pwMin[p] = min_v;
                        pwWidth[p] = max_v - min_v;
                        totalLen += (max_v - min_v) * svpar.pieceLength;
                    }
                    
                    Batch_SV_Data[k][v_idx].length = totalLen;
                    Batch_SV_Data[k][v_idx].pieceWiseMin = pwMin;
                    Batch_SV_Data[k][v_idx].pieceWiseWidth = pwWidth;
                    
                    unsigned char *padded = (unsigned char*)calloc(totalLen, sizeof(unsigned char));
                    unsigned char *transposed = (unsigned char*)malloc(totalLen * sizeof(unsigned char));
                    
                    unsigned char *ptr_raw = AVal_arr[r_y][r_x].val;
                    unsigned char *ptr_pad = padded;
                    
                    for(int p=0; p<NViews; p++) {
                        int set_idx = p / svpar.pieceLength;
                        int n_pad_pre = (int)ACol_arr[r_y][r_x].minIndex[p] - (int)pwMin[set_idx] - (int)bandMin[p];
                        ptr_pad += n_pad_pre; 
                        
                        int count = ACol_arr[r_y][r_x].countTheta[p];
                        memcpy(ptr_pad, ptr_raw, count);
                        ptr_pad += count;
                        ptr_raw += count;
                        
                        int n_pad_post = (int)pwMin[set_idx] + (int)pwWidth[set_idx] - (int)ACol_arr[r_y][r_x].minIndex[p] - count + (int)bandMin[p];
                        ptr_pad += n_pad_post;
                    }
                    
                    unsigned char *ptr_src = padded;
                    unsigned char *ptr_dst = transposed;
                    for(int p=0; p<NViewSets; p++) {
                        int w = pwWidth[p];
                        int h = svpar.pieceLength;
                        for(int col=0; col<w; col++) {
                            for(int row=0; row<h; row++) {
                                ptr_dst[col*h + row] = ptr_src[row*w + col];
                            }
                        }
                        ptr_src += w*h;
                        ptr_dst += w*h;
                    }
                    
                    free(padded);
                    Batch_SV_Data[k][v_idx].val = transposed;
                }
            }
        } // Parallel region end

        // Write to file
        int SV_Vol_Count = (2*SVLength+1)*(2*SVLength+1);
        for(int k=0; k<current_batch_count; k++) {
            int sv_idx = batch_sv_indices[k];

            fwrite(svpar.bandMinMap[sv_idx].bandMin, sizeof(channel_t), NViews, fp);
            fwrite(svpar.bandMaxMap[sv_idx].bandMax, sizeof(channel_t), NViews, fp);

            for(int v=0; v<SV_Vol_Count; v++) {
                int len = Batch_SV_Data[k][v].length;
                fwrite(&len, sizeof(int), 1, fp);
                if(len > 0) {
                    fwrite(Batch_SV_Data[k][v].val, sizeof(unsigned char), len, fp);
                    fwrite(Batch_SV_Data[k][v].pieceWiseMin, sizeof(channel_t), NViewSets, fp);
                    fwrite(Batch_SV_Data[k][v].pieceWiseWidth, sizeof(channel_t), NViewSets, fp);
                    
                    free(Batch_SV_Data[k][v].val);
                    free(Batch_SV_Data[k][v].pieceWiseMin);
                    free(Batch_SV_Data[k][v].pieceWiseWidth);
                }
            }
            free(Batch_SV_Data[k]);
        }
        free(Batch_SV_Data);

        for(i = y_min; i < y_max; i++) {
            for(j = 0; j < Nx; j++) {
                if(ACol_arr[i][j].n_index > 0) {
                    free(ACol_arr[i][j].countTheta);
                    free(ACol_arr[i][j].minIndex);
                    free(AVal_arr[i][j].val);
                }
            }
            free(ACol_arr[i]);
            free(AVal_arr[i]);
        }
        free(ACol_arr);
        free(AVal_arr);
        free(batch_sv_indices);
    } 

    if(verboseLevel) fprintf(stdout, "Finalizing file (Aval_max)...\n");
    fseek(fp, 0, SEEK_END);
    fwrite(Aval_max_ptr, sizeof(float), Nx*Ny, fp);

    fclose(fp);
    
    for(i=0; i<Nsv; i++) {
        free(svpar.bandMinMap[i].bandMin);
        free(svpar.bandMaxMap[i].bandMax);
    }
    free(svpar.bandMinMap);
    free(svpar.bandMaxMap);
    free(Aval_max_ptr);
    free(order);
    free_img((void **)pix_prof);
    
    if(verboseLevel) fprintf(stdout, "Streaming computation complete.\n");
}
