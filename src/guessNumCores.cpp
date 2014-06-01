#include "config.hpp"
#include "guessNumCores.hpp"

#include <cstddef> // size_t

#include <external/io.h>

using std::uint32_t;
using std::size_t;



// windows adapted from
//   http://scriptionary.com/2010/09/01/counting-processor-cores-and-threads/
//   http://msdn.microsoft.com/en-us/library/ms683194%28v=VS.85%29.aspx
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>

#  ifndef _W64
#    include "glpi.h"
#  endif

typedef BOOL (WINAPI* glpiFunction)(PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, PDWORD);

namespace cibart {
  void guessNumCores(uint32_t* numPhyiscalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhyiscalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    glpiFunction getLogicalProcessorInformation = (glpiFunction) GetProcAddress(GetModuleHandle(TEXT("kernel32")), "GetLogicalProcessorInformation");
    if (NULL == glpi) {
      issueWarning("GetLogicalProcessorInformation is not supported on this OS.");
      return;
    }
    
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = NULL;
    DWORD bufferLength = 0;
    
    // retrieve appropriate buffer length
    getLogicalProcessorInformation(buffer, &bufferLength);
    if (GetLastError() != ERROR_INSUFFICIENT_BUFFER) {
      issueWarning("Unable to determine buffer length for GetLogicalProcessorInformation: %d", GetLastError());
      return;
    }
    
    buffer = new SYSTEM_LOGICAL_PROCESSOR_INFORMATION[bufferLength];
    if (getLogicalProcessorInformation(buffer, &bufferLength) == false) {
      delete [] buffer;
      issueWarning("Unable to allocate buffer length for GetLogicalProcessorInformation");
      return;
    }
    
    size_t numElements = bufferLength / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
    
    for (size_t i = 0; i < numElements; ++i) {
      if (buffer[i].Relationship != RelationProcessorCore) continue;
      
      *numPhyiscalProcessors += 1;
      
      if (buffer[i].ProcessorCore.Flags == 1) {
        int mask = 0x1;
        size_t numBitsInInt = 8 * sizeof(ULONG_PTR);
        for (size_t j = 0; j < numBitsInInt; ++j) {
          if (buffer[i].ProcessorMask & mask == true) *numLogicalProcessors += 1;
          mask <<= 1;
        }
      }
    }
    
    delete [] buffer;
  }
}

#elif defined(__MACH__)
#  include <sys/types.h>
#  include <sys/sysctl.h>

namespace cibart {
  void guessNumCores(uint32_t* numPhyiscalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhyiscalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    int query[2];
    size_t queryLength = 2;
    
    size_t bufferLength = sizeof(uint32_t);
    
    // these have apparently been around since 2003 so I'm not sure when/if they could fail
    if (sysctlnametomib("hw.physicalcpu", query, &queryLength) != -1) {
      sysctl(query, queryLength, numPhyiscalProcessorsPtr, &bufferLength, NULL, 0);
    }
    
    queryLength = 2;
    bufferLength = sizeof(uint32_t);
    
    if (sysctlnametomib("hw.logicalcpu", query, &queryLength) != -1) {
      sysctl(query, queryLength, numLogicalProcessorsPtr, &bufferLength, NULL, 0);
    }
    
    bufferLength = sizeof(uint32_t);
    
    if (*numLogicalProcessorsPtr == 0) {
      query[0] = CTL_HW;
#  ifdef HW_AVAILCPU
      query[1] = HW_AVAILCPU;
      if (sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0) == -1 || *numLogicalProcessorsPtr == 0) {
#  endif
        query[1] = HW_NCPU;
        sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0);
      }
#  ifdef HW_AVAILCPU
    }
#  endif
  }
}

// the only way to reliably determine the number of physical cores
// is to parse /proc/cpuinfo; sysconf(_SC_NPROCESSORS_CONF) returns
// hyper as well
#elif defined(__linux__)
#  include <fcntl.h>    // open
#  include <unistd.h>   // close, lseek
#  include <sys/mman.h> // mmap
#  include <string.h>   // srtncmp
#  include <errno.h>

#  include <external/alloca.h>

#  include <map>
#  include <vector>
#  include <utility>

namespace {
  typedef std::map<uint32_t, uint32_t> CoreMap;
  typedef std::pair<uint32_t, uint32_t> CorePair;
  struct Processor {
    uint32_t physicalId;
    CoreMap coreIdThreadMap;
  };
  typedef std::map<uint32_t, Processor*> CPUMap;
  typedef std::pair<uint32_t, Processor*> CPUPair;
  
  enum LineType {
    PROCESSOR = 0,
    PHYSICAL_ID,
    CORE_ID,
    NA
  };
  
#  define INVALID_ID ((uint32_t) -1)
  
  void parseRecord(CPUMap& processorIdMap, uint32_t& currentPhysicalId, uint32_t& currentCoreId, uint32_t& currentProcessorId, uint32_t newId)
  {
    Processor* currentProcessor;
    
    // new records are introduced by processor id lines
    if (currentProcessorId != newId && currentProcessorId != INVALID_ID) {
      
      if (currentPhysicalId == INVALID_ID) {
        // only should happen when there is a single processor with one core
        // and no hyperthreading
        currentProcessor = new Processor;
        currentProcessor->physicalId = 0;
        currentProcessor->coreIdThreadMap.insert(CorePair(0, 1));
        
        processorIdMap.insert(CPUPair(currentProcessorId, currentProcessor));
      } else {
        // this should never happen
        if (currentCoreId == INVALID_ID) currentCoreId = 0;
        
        CPUMap::iterator processorLookup = processorIdMap.find(currentPhysicalId);
        if (processorLookup == processorIdMap.end()) {
          currentProcessor = new Processor;
          currentProcessor->physicalId = currentPhysicalId;
          currentProcessor->coreIdThreadMap.insert(CorePair(currentCoreId, 1));
          processorIdMap.insert(CPUPair(currentPhysicalId, currentProcessor));
        } else {
          currentProcessor = processorLookup->second;
          CoreMap::iterator coreLookup = currentProcessor->coreIdThreadMap.find(currentCoreId);
          if (coreLookup == currentProcessor->coreIdThreadMap.end()) {
            currentProcessor->coreIdThreadMap.insert(CorePair(currentCoreId, 1));
          } else {
            coreLookup->second += 1;
          }
        }
      }
      
      
      currentProcessorId = newId;
      currentPhysicalId = INVALID_ID;
      currentCoreId = INVALID_ID;
    } else {
      // first line, usually
      currentProcessorId = newId;
    }
  }

#  define ERROR_BUFFER_LENGTH 1024
  bool parseProcCPUInfo(std::vector < Processor* >& result)
  {
    char errorBuffer[ERROR_BUFFER_LENGTH];
    
    int fd = open("/proc/cpuinfo", O_RDONLY);
    if (fd == -1) {
      strerror_r(errno, errorBuffer, ERROR_BUFFER_LENGTH);
      issueWarning("unable to open /proc/cpuinfo: %s\n", errorBuffer);
      return false;
    }
    
    off_t fileLength = lseek(fd, 0, SEEK_END);
    if (fileLength == ((off_t) - 1)) {
      strerror_r(errno, errorBuffer, ERROR_BUFFER_LENGTH);
      issueWarning("unable to seek in /proc/cpuinfo: %s\n", errorBuffer);
      close(fd);
      return false;
    }
    lseek(fd, 0, SEEK_SET);
    
    void* v_cpuInfo = mmap(0, fileLength, PROT_READ, MAP_PRIVATE, fd, 0);
    if (v_cpuInfo == MAP_FAILED) {
      strerror_r(errno, errorBuffer, ERROR_BUFFER_LENGTH);
      issueWarning("unable to map /proc/cpuinfo: %s\n", errorBuffer);
      close(fd);
      return false;
    }
    
    char* cpuInfo = (char*) v_cpuInfo;
    
    off_t offset = 0;
    
    uint32_t currentProcessorId = INVALID_ID;
    uint32_t currentPhysicalId  = INVALID_ID;
    uint32_t currentCoreId      = INVALID_ID;
    CPUMap processorIdMap;
    LineType lineType;
    
    while (offset < fileLength) {
      lineType = NA;
      
      if ((cpuInfo[offset] == 'p' || cpuInfo[offset] == 'P') && fileLength - offset >= 9 &&
          strncmp(cpuInfo + offset + 1, "rocessor", 8) == 0) {
        lineType = PROCESSOR;
        offset += 9;
      } else if ((cpuInfo[offset] == 'p' || cpuInfo[offset] == 'P') && fileLength - offset >= 11 &&
                  strncmp(cpuInfo + offset + 1, "hysical id", 10) == 0) {
        lineType = PHYSICAL_ID;
        offset += 11;
      } else if ((cpuInfo[offset] == 'c' || cpuInfo[offset] == 'C') && fileLength - offset >= 7 &&
                  strncmp(cpuInfo + offset + 1, "ore id", 6) == 0) {
        lineType = CORE_ID;
        offset += 7;
      }
      
      if (lineType == NA) {
        while (offset < fileLength && cpuInfo[offset] != '\n') ++offset;
        ++offset;
        if (offset >= fileLength) break;
        continue;
      }
      
      while (offset < fileLength && cpuInfo[offset] != ':' && cpuInfo[offset] != '\n') ++offset;
      if (cpuInfo[offset] == '\n') {
        ++offset; continue;
      }
      if (offset == fileLength) break;
      ++offset;
      
      while (offset < fileLength && cpuInfo[offset] < '0' && cpuInfo[offset] > '9' && cpuInfo[offset] != '\n') ++offset;
      if (cpuInfo[offset] == '\n') {
        ++offset; continue;
      }
      if (offset == fileLength) break;
      
      off_t endOfNumber = offset + 1;
      while (endOfNumber < fileLength && cpuInfo[endOfNumber] >= '0' && cpuInfo[endOfNumber] <= '9') ++endOfNumber;
      
      char* buffer = ext_stackAllocate(endOfNumber - offset + 1, char);
      memcpy(buffer, (const char*) cpuInfo + offset, endOfNumber - offset);
      buffer[endOfNumber - offset] = '\0';
      
      long parsedInt = strtol(buffer, NULL, 10);
      ext_stackFree(buffer);
      
      offset = endOfNumber;
      
      
      switch (lineType) {
        case PROCESSOR:
          parseRecord(processorIdMap, currentPhysicalId, currentCoreId, currentProcessorId, parsedInt);
          break;
        case PHYSICAL_ID:
          currentPhysicalId = parsedInt;
          break;
        case CORE_ID:
          currentCoreId = parsedInt;
          break;
        default:
          break;
      }
      
      while (offset < fileLength && cpuInfo[offset] != '\n') ++offset;
      ++offset;
      if (offset >= fileLength) break;
    }
    
    parseRecord(processorIdMap, currentPhysicalId, currentCoreId, currentProcessorId, INVALID_ID);
    
    for (CPUMap::iterator it = processorIdMap.begin(); it != processorIdMap.end(); ++it) {
      result.push_back(it->second);
    }
    
    munmap(v_cpuInfo, fileLength);
    close(fd);
    
    return true;
  }
}

namespace cibart {
  void guessNumCores(uint32_t* numPhyiscalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhyiscalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    std::vector<Processor*> cpuInfo;
    if (parseProcCPUInfo(cpuInfo) == true) {
      for (size_t i = 0; i < cpuInfo.size(); ++i) {
        Processor* processor = cpuInfo[i];
        *numPhyiscalProcessorsPtr += processor->coreIdThreadMap.size();
        for (CoreMap::iterator it = processor->coreIdThreadMap.begin(); it != processor->coreIdThreadMap.end(); ++it) {
          *numLogicalProcessorsPtr += it->second;
        }
      }
    }
#  if defined(_SC_NPROCESSORS_ONLN)
    else {
      *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_ONLN);
      if (*numLogicalProcessorsPtr < 1) *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_CONF);
    }
#  endif
    
    for (size_t i = 0; i < cpuInfo.size(); ++i) {
      Processor* processor = cpuInfo[i];
      delete processor;
    }
  }
}

#else // freeball it

#  include <unistd.h>
#  ifdef _SC_NPROCESSORS_ONLN
#    define USE_SYSCONF
#  elif defined(HAVE_SYS_SYSCTL_H)
#    include <sys/param.h>
#    include <sys/sysctl.h>
#    define USE_SYSCTL
#  endif


namespace cibart {
  void guessNumCores(uint32_t* numPhyiscalProcessorsPtr, uint32_t* numLogicalProcessorsPtr) {
    *numPhyiscalProcessorsPtr = 0;
    *numLogicalProcessorsPtr  = 0;
    
    
#  if defined(USE_SYSCONF)
    *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_ONLN);
    if (*numLogicalProcessorsPtr < 1) {
#    ifdef _SC_NPROCESSORS_CONF
      *numLogicalProcessorsPtr = sysconf(_SC_NPROCESSORS_CONF);
      if (result < 1) result = -1;
#    else
      result = -1;
#    endif
    }
    
#  elif defined(USE_SYSCTL)
    int query[2];
    size_t bufferLength = sizeof(uint32_t);
    
    query[0] = CTL_HW;
#    ifdef HW_AVAILCPU
    query[1] = HW_AVAILCPU;
    if (sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0) == -1 ||  numCores < 1) {
#    endif
      key[1] = HW_NCPU;
      sysctl(query, 2, numLogicalProcessorsPtr, &bufferLength, NULL, 0);
#    ifdef HW_AVAILCPU
    }
#    endif
#  endif
  }
}

#endif
