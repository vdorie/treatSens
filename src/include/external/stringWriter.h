#ifndef EXTERNAL_STRING_WRITER_H
#define EXTERNAL_STRING_WRITER_H

// wraps a character buffer and lengthens it as needed

#include "stddef.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EXT_SWR_DEFAULT_BUFFER_SIZE 1024

typedef struct ext_stringWriter
{
  char* buffer;
  ext_size_t length;
  ext_size_t pos;
} ext_stringWriter;

ext_stringWriter* ext_swr_create(ext_size_t bufferLength);
void ext_swr_destroy(ext_stringWriter* writer);
int ext_swr_initialize(ext_stringWriter* writer, ext_size_t bufferLength);
void ext_swr_invalidate(ext_stringWriter* writer);
int ext_swr_reset(ext_stringWriter* writer);

int ext_swr_writeChar(ext_stringWriter* writer, char c);
int ext_swr_writeString(ext_stringWriter* writer, const char* c, ext_size_t length);
int ext_swr_write32BitInteger(ext_stringWriter* writer, int32_t i);

#ifdef __cplusplus
}
#endif

#endif // EXTERNAL_STRING_WRITER_H
