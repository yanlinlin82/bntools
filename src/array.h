#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <stdlib.h>
#include <errno.h>

#define array(type) struct { size_t capacity, size; type *data; }
typedef array(void) array_t;

#define array_init(var) \
	do { (var).capacity = (var).size = 0; (var).data = NULL; } while (0)

#define array_free(var) \
	do { free((var).data); (var).capacity = (var).size = 0; (var).data = NULL; } while (0)

static inline size_t array_new_capacity_(size_t capacity, size_t size)
{
	const size_t MIN_INCREMENT = 16;
	const size_t MAX_INCREMENT = 1024;
	while (capacity < size) {
		if (capacity < MIN_INCREMENT) {
			capacity += MIN_INCREMENT;
		} else if (capacity > MAX_INCREMENT) {
			capacity += MAX_INCREMENT;
		} else {
			capacity += capacity;
		}
	}
	return capacity;
}

#define array_reserve(var, new_size) \
	({ \
		int res_ = 0; \
		size_t size_ = (new_size); \
		if (size_ > (var).capacity) { \
			size_t capacity_ = array_new_capacity_((var).capacity, size_); \
			typeof((var).data) p_ = malloc(sizeof(*p_) * capacity_); \
			if (!p_) { \
				res_ = -ENOMEM; \
			} else { \
				if ((var).size > 0) { \
					memcpy(p_, (var).data, sizeof(*p_) * (var).size); \
				} \
				memset(p_ + (var).size, 0, sizeof(*p_) * (capacity_ - (var).size)); \
				free((var).data); \
				(var).data = p_; \
				(var).capacity = capacity_; \
			} \
		} \
		res_; \
	})

#endif /* __ARRAY_H__ */
