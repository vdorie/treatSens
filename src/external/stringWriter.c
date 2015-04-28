#include <external/stringWriter.h>

#include <errno.h>
#include <stdlib.h> // malloc
#include <string.h> // memcpy
#include <stdio.h>  // snprintf

#define INT_BUFFER_LENGTH 16

static int reallocate(ext_stringWriter* writer) {
  char* temp = (char*) malloc(2 * writer->length);
  if (temp == NULL) return errno;
  
  memcpy(temp, (const char*) writer->buffer, writer->length * sizeof(char));
  writer->length *= 2;
  free(writer->buffer);
  writer->buffer = temp;
  
  return 0;
}

ext_stringWriter* ext_swr_create(size_t bufferLength)
{
  ext_stringWriter* result = (ext_stringWriter*) malloc(sizeof(ext_stringWriter));
  if (result == NULL) return NULL;
  
  int errorCode = ext_swr_initialize(result, bufferLength);
  if (errorCode != 0) {
    free(result);
    return NULL;
  }
  
  return result;
}

void ext_swr_destroy(ext_stringWriter* writer)
{
  if (writer != NULL) {
    ext_swr_invalidate(writer);
    free(writer);
  }
}

int ext_swr_initialize(ext_stringWriter* writer, size_t bufferLength)
{
  if (writer == NULL) return EINVAL;
  
  writer->buffer = (char*) malloc(bufferLength * sizeof(char));
  if (writer->buffer == NULL) return errno;
  
  writer->length = bufferLength;
  writer->pos = 0;
  
  return 0;
}

int ext_swr_reset(ext_stringWriter* writer) 
{
  if (writer == NULL) return EINVAL;
  
  writer->pos = 0;
  
  return 0;
}

void ext_swr_invalidate(ext_stringWriter* writer)
{
  if (writer != NULL && writer->buffer != NULL) free(writer->buffer);
}

int ext_swr_writeChar(ext_stringWriter* writer, char c)
{
  if (writer == NULL) return EINVAL;
  
  writer->buffer[writer->pos++] = c;
  
  if (writer->pos >= writer->length) return reallocate(writer);
  
  return 0;
}

int ext_swr_writeString(ext_stringWriter* writer, const char* s, size_t length)
{
  if (writer == NULL) return EINVAL;
  
  int errorCode = 0;
  
  if (writer->pos + length >= writer->length && (errorCode = reallocate(writer)) != 0) return errorCode;
  
  memcpy(writer->buffer + writer->pos, s, length * sizeof(char));
  writer->pos += length;
  
  return 0;
}

int ext_swr_write32BitInteger(ext_stringWriter* writer, int32_t i)
{
  if (writer == NULL) return EINVAL;
  
  char intBuffer[INT_BUFFER_LENGTH];
  
  // 2^32 == 4,294,967,296 <==> 10 chars, 11 with sign; can't imagine how this can ever fail
  
  int bytesWritten = snprintf(intBuffer, INT_BUFFER_LENGTH, "%d", i);

  return ext_swr_writeString(writer, intBuffer, (size_t) bytesWritten);
}
