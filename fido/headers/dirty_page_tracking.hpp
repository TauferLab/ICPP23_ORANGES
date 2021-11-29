#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <chrono>
#include <inttypes.h>

typedef struct {
  uint64_t pfn : 55;
  uint64_t soft_dirty : 1;
  uint64_t exclusive_mapped : 1;
  uint64_t zeros : 4;
  uint64_t file_page : 1;
  uint64_t swapped : 1;
  uint64_t present : 1;
} PagemapEntry;

/* Parse the pagemap entry for the given virtual address.
 *
 * @param[out] entry      the parsed entry
 * @param[in]  pagemap_fd file descriptor to an open /proc/pid/pagemap file
 * @param[in]  vaddr      virtual address to get entry for
 * @return 0 for success, 1 for failure
 */
int pagemap_get_entry(PagemapEntry& entry, int pagemap_fd, uintptr_t vaddr) {
  size_t nread;
  ssize_t ret;
  uint64_t data;
  size_t page_size = sysconf(_SC_PAGE_SIZE);
  
  nread = 0;
  while(nread < sizeof(data)) {
    size_t count = sizeof(data) - nread;
    off_t offset = (vaddr / page_size) * sizeof(data) + nread;
    ret = pread(pagemap_fd, reinterpret_cast<uint8_t*>(&data) + nread, count, offset);
    nread += ret;
    if(ret <=0) {
      return -1;
    }
  }
  entry.pfn = data & ((static_cast<uint64_t>(1) << 55) - 1);
  entry.soft_dirty = (data >> 55) & 1;
  entry.exclusive_mapped = (data >> 56) & 1;
  entry.file_page = (data >> 61) & 1;
  entry.swapped = (data >> 62) & 1;
  entry.present = (data >> 63) & 1;
  return 0;
}

bool address_dirty(pid_t pid, uintptr_t vaddr) {
  char pagemap_file[BUFSIZ];
  int pagemap_fd;
  
  snprintf(pagemap_file, sizeof(pagemap_file), "/proc/%ju/pagemap", static_cast<uintmax_t>(pid));
  pagemap_fd = open(pagemap_file, O_RDONLY);
  if(pagemap_fd < 0) {
    return 1;
  }
  PagemapEntry entry;
  if(pagemap_get_entry(entry, pagemap_fd, vaddr)) {
    return 1;
  }
  close(pagemap_fd);
  printf("0x%016" PRIXPTR "\n", (int64_t)(static_cast<int64_t>(1) << 55));
  printf("0x%016" PRIXPTR "\n", entry);
  if(entry.soft_dirty)
    return true;
  return false;
}

bool address_range_dirty(pid_t pid, uintptr_t start_addr, uintptr_t end_addr) {
  char pagemap_file[BUFSIZ];
  int pagemap_fd;
  uintptr_t current = start_addr;
  uint64_t page_size = sysconf(_SC_PAGE_SIZE);

  snprintf(pagemap_file, sizeof(pagemap_file), "/proc/%ju/pagemap", static_cast<uintmax_t>(pid));
  pagemap_fd = open(pagemap_file, O_RDONLY);
  if(pagemap_fd < 0) {
    printf("Failed to open pagemap file\n");
    return true;
  }

  int page_counter = 0;
  PagemapEntry entry;
  while(current < end_addr) {
    int ret = pagemap_get_entry(entry, pagemap_fd, current);
    if(ret)
      printf("Something bad happened when getting the pagemap entry\n");
    if(entry.soft_dirty) {
      return true;
    }
    printf("0x%016" PRIXPTR "\n", entry);
    current += page_size;
  }
  int ret = pagemap_get_entry(entry, pagemap_fd, end_addr);
  if(ret)
    printf("Something bad happened when getting the pagemap entry\n");
  if(entry.soft_dirty) {
    printf("0x%016" PRIXPTR "\n", entry);
    return true;
  }

  close(pagemap_fd);
  return false;
}

bool get_dirty_pages(pid_t pid, uintptr_t start_addr, uintptr_t end_addr, std::vector<uint64_t>& page_list) {
  char pagemap_file[BUFSIZ];
  int pagemap_fd;
  uintptr_t current = start_addr;
  uint64_t page_size = sysconf(_SC_PAGE_SIZE);

  page_list.clear();

  snprintf(pagemap_file, sizeof(pagemap_file), "/proc/%ju/pagemap", static_cast<uintmax_t>(pid));
  pagemap_fd = open(pagemap_file, O_RDONLY);
  if(pagemap_fd < 0) {
    printf("Failed to open pagemap file\n");
    return true;
  }

  int page_counter = 0;
  PagemapEntry entry;
  while(current < end_addr) {
    int ret = pagemap_get_entry(entry, pagemap_fd, current);
    if(ret)
      printf("Something bad happened when getting the pagemap entry\n");
    if(entry.soft_dirty) {
//      page_list[page_counter++] = (current/page_size);
      page_list.push_back(current/page_size);
      page_counter += 1;
    }
//    printf("0x%016" PRIXPTR "\n", entry);
    current += page_size;
  }
  int ret = pagemap_get_entry(entry, pagemap_fd, end_addr);
  if(ret)
    printf("Something bad happened when getting the pagemap entry\n");
  if(entry.soft_dirty) {
//    page_list[page_counter] = (end_addr/page_size);
    page_list.push_back(end_addr/page_size);
    page_counter += 1;
//    printf("0x%016" PRIXPTR "\n", entry);
  }

  close(pagemap_fd);
  if(page_counter > 0) {
    return true;
  } else {
    return false;
  }
}

void reset_dirty_bit(pid_t pid) {
  int fd, ret;
  char clear_command[] = "4";
  char clear_file[BUFSIZ];
  snprintf(clear_file, sizeof(clear_file), "/proc/%ju/clear_refs", static_cast<uintmax_t>(pid));
  fd = open(clear_file, O_WRONLY);
  ret = write(fd, clear_command, sizeof(clear_command));
  close(fd);
  return;
}

