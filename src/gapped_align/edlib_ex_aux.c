#include "edlib_ex_aux.h"

EdlibAlignMatrix*
new_EdlibAlignMatrix()
{
    EdlibAlignMatrix* data = (EdlibAlignMatrix*)malloc(sizeof(EdlibAlignMatrix));
    data->Ps = (Word*)malloc(sizeof(Word) * MaxNumBlocks * MaxSeqSize);
    data->Ms = (Word*)malloc(sizeof(Word) * MaxNumBlocks * MaxSeqSize);
    data->scores = (int*)malloc(sizeof(int) * MaxNumBlocks * MaxSeqSize);
    data->first_blocks = (int*)malloc(sizeof(int) * MaxSeqSize);
    data->last_blocks = (int*)malloc(sizeof(int) * MaxSeqSize);

    return data;
}

EdlibAlignMatrix*
free_EdlibAlignMatrix(EdlibAlignMatrix* data)
{
    free(data->Ps);
    free(data->Ms);
    free(data->scores);
    free(data->first_blocks);
    free(data->last_blocks);
    free(data);
    return 0;
}

EdlibAlignResult*
new_EdlibAlignResult()
{
    EdlibAlignResult* result = (EdlibAlignResult*)malloc(sizeof(EdlibAlignResult));
    kv_init(result->end_locations);
    kv_init(result->start_locations);
    kv_init(result->alignment);
    return result;
}

EdlibAlignResult*
free_EdlibAlignResult(EdlibAlignResult* result)
{
    kv_destroy(result->end_locations);
    kv_destroy(result->start_locations);
    kv_destroy(result->alignment);
    free(result);
    return 0;
}

void
clear_EdlibAlignResult(EdlibAlignResult* result)
{
    result->edit_distance = -1;
    kv_clear(result->end_locations);
    kv_clear(result->start_locations);
    kv_clear(result->alignment);
}

EdlibAlignData*
new_EdlibAlignData()
{
    EdlibAlignData* data = (EdlibAlignData*)malloc(sizeof(EdlibAlignData));
    data->align_matrix = new_EdlibAlignMatrix();
    data->result = new_EdlibAlignResult();
    data->blocks = (EdlibBlock*)malloc(sizeof(EdlibBlock) * MaxNumBlocks);
	data->peq = (Word*)malloc(sizeof(Word) * (AlphabetSize + 1) * MaxNumBlocks);
    kv_init(data->cigar);
    return data;
}

EdlibAlignData*
free_EdlibAlignData(EdlibAlignData* data)
{
    data->align_matrix = free_EdlibAlignMatrix(data->align_matrix);
    data->result = free_EdlibAlignResult(data->result);
    free(data->blocks);
	free(data->peq);
    kv_destroy(data->cigar);
    free(data);
    return 0;
}
