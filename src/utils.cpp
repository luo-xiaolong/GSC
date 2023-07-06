// *******************************************************************************************
// This file is a part of VCFShark software distributed under GNU GPL 3 licence.
// The homepage of the VCFShark project is https://github.com/refresh-bio/VCFShark
//
// Authors: Sebastian Deorowicz, Agnieszka Danek, Marek Kokot
// Version: 1.1
// Date   : 2021-02-18
// *******************************************************************************************

#include "utils.h"
#include <iostream>
#include <memory>
#include <sstream>
#include <algorithm>

#ifdef OUR_STRTOL
// *****************************************************************************************
#ifdef __APPLE__
long int strtol(const char* str, char** endptr, int base)
#else
long int strtol(const char* str, char** endptr, int base) noexcept
#endif
{
	if (base != 10)
	{
		std::cerr << "unsuported base " << base << std::endl;
		fflush(stdout);
		exit(1);
	}

	long int val = 0;
	char *p = (char*)str;
	bool is_negative = false;

	if (*p == '-')
	{
		is_negative = true;
		++p;
	}

	while (*p >= '0' && *p <= '9') 
	{
		val = val * 10 + (*p++ - '0');
	}

	if (endptr)
		*endptr = p;

	return is_negative ? -val : val;
}
#endif
// ************************************************************************************
void append_str(vector<uint8_t>& v_comp, const string& x) {
    size_t old_size = v_comp.size();
    size_t new_size = old_size + x.size() + 1; // 加 1 表示空字符
    v_comp.reserve(new_size);
    copy(x.begin(), x.end(), std::back_inserter(v_comp));
    v_comp.emplace_back('\0'); // 在末尾添加空字符

}
void read_str(const vector<uint8_t>& v_comp, size_t& pos, string& x)
{
	x.clear();
    auto null_pos = find(v_comp.begin() + pos, v_comp.end(),'\0');
    x.assign(v_comp.begin() + pos, null_pos);
    pos = null_pos - v_comp.begin() + 1;

}
void my_merge(vector<sblock> &a, int start, int mid, int end)
{
	int size_left = mid - start + 1; // 左半段数量
	int size_right = end - mid;		 // 右半段数量
	vector<sblock> left;			 // 加1是为了多加一个哨兵位
	vector<sblock> right;
	int i = 0, j = 0, k = 0; // 循环计数变量
	for (i = 0; i < size_left; ++i)
	{
		left.emplace_back(a[start + i]);
	}
	for (i = 0; i < size_right; ++i)
	{
		right.emplace_back(a[mid + 1 + i]);
	}
	// 加哨兵
	sblock max(MAX, 0, "shaobing");
	left.emplace_back(max);
	right.emplace_back(max);
	// 按序生成新的数组
	for (k = start, i = 0, j = 0; k <= end; ++k)
	{
		if (left[i].val < right[j].val)
		{
			a[k] = left[i];
			++i;
		}
		else if (left[i].val == right[j].val)
		{
			if (atoi(left[i].s_ref.c_str()) <= atoi(right[j].s_ref.c_str()))
			{
				a[k] = left[i];
				++i;
			}
			else
			{
				a[k] = right[j];
				++j;
			}
		}
		else
		{
			a[k] = right[j];
			++j;
		}
	}
	left.clear();
	right.clear();
}
// *******************************************************************************************************************************
void my_merge_sort_r(vector<sblock> &a, int start, int end)
{
	if (start >= end)
	{
		return;
	}
	int mid = start + ((end - start) >> 1); // 中间元素下标
											// 两段数据分开进行排序
	my_merge_sort_r(a, start, mid);
	my_merge_sort_r(a, mid + 1, end);
	// 合并
	my_merge(a, start, mid, end);
}
// *******************************************************************************************************************************
void my_merge_sort(vector<sblock> &a, int len)
{
	if (a.empty() || len <= 1)
	{
		return;
	}
	// 递归进行归并排序
	my_merge_sort_r(a, 0, len - 1);
}

// // *******************************************************************************************************************************
// void insertion_sort(vector<sblock>& arr) {
//     int n = arr.size();
//     for (int i = 1; i < n; ++i) {
//         sblock key = arr.at(i);
//         int j = i - 1;
//         while (j >= 0 && (arr.at(j).val > key.val || (arr.at(j).val == key.val && atoi(arr.at(j).s_ref.c_str()) > atoi(key.s_ref.c_str())))) {
//             arr.at(j + 1) = arr.at(j);
//             j = j - 1;
//         }
//         arr.at(j + 1) = key;
//     }
// }

// void merge(vector<sblock>& arr, int l, int m, int r) {
//     int n1 = m - l + 1;
//     int n2 = r - m;
 
//     vector<sblock> L, R;
 
//     for (int i = 0; i < n1; i++) {
//         L.emplace_back(arr.at(l + i)) ;
//     }
//     for (int j = 0; j < n2; j++) {
//         R.emplace_back(arr.at(m + 1 + j));
//     }
 
//     int i = 0, j = 0, k = l;
 
//     while (i < n1 && j < n2) {
//         if (L.at(i).val < R.at(j).val || (L.at(i).val == R.at(j).val && atoi(L.at(i).s_ref.c_str()) <= atoi(R.at(j).s_ref.c_str()))) {
//             arr.at(k) = L.at(i);
//             i++;
//         }
//         else {
//             arr.at(k) = R.at(j);
//             j++;
//         }
//         k++;
//     }
 
//     while (i < n1) {
//         arr.at(k) = L.at(i);
//         i++;
//         k++;
//     }
 
//     while (j < n2) {
//         arr.at(k) = R.at(j);
//         j++;
//         k++;
//     }
// }

// void merge_sort(vector<sblock>& arr, int l, int r) {
//     if (l >= r) {
//         return;
//     }
//     if (r - l + 1 <= 5) { // 对于小于等于5个元素的子数组，采用插入排序
//         insertion_sort(arr);
//         return;
//     }
//     int m = l + (r - l) / 2;
//     merge_sort(arr, l, m);
//     merge_sort(arr, m + 1, r);
//     merge(arr, l, m, r);
// }

// void my_merge_sort(vector<sblock>& arr) {
//     int n = arr.size();
//     merge_sort(arr, 0, n - 1);
// }
static inline uint8_t *bcf_unpack_fmt_core1(uint8_t *ptr, int n_sample, bcf_fmt_t *fmt)
{

    uint8_t *ptr_start = ptr;
    fmt->id = bcf_dec_typed_int1(ptr, &ptr);
    fmt->n = bcf_dec_size(ptr, &ptr, &fmt->type);
    fmt->size = fmt->n << bcf_type_shift[fmt->type];
    fmt->p = ptr;
    fmt->p_off  = ptr - ptr_start;
    fmt->p_free = 0;
    ptr += n_sample * fmt->size;
    fmt->p_len = ptr - fmt->p;
    return ptr;
}

// Based on function from htslib: bcf_update_format((hdr),(line),"GT",(gts),(n),BCF_HT_INT)
int bcf_update_genotypes_fast(const bcf_hdr_t *hdr, bcf1_t *line, kstring_t &str)
{
	bcf_fmt_t *fmt = NULL;

	line->n_sample = bcf_hdr_nsamples(hdr);

	// bcf_enc_vint(&str, n, (int32_t*)values, nps); 
	// bcf_unpack(line, BCF_UN_ALL);

    line->n_fmt++;
	// cout<<line->n_fmt<<endl;
	// 	int i;
	// if(!a)
	// 	for (i=0; i<line->n_fmt; i++){
	// 	cout<<line->d.fmt[i].id<<endl;
	// 	// if ( line->d.fmt[i].id == (int)str.s[1] ) break;

	// }
	// a=1;
    hts_expand00(bcf_fmt_t, line->n_fmt, line->d.m_fmt, line->d.fmt);
	 
    fmt = &line->d.fmt[line->n_fmt-1];
	
    bcf_unpack_fmt_core1((uint8_t*)str.s, line->n_sample, fmt);
	// cout<<fmt->id<<endl;
    line->d.indiv_dirty = 1;
    fmt->p_free = 0;
    line->unpacked |= BCF_UN_FMT;

    return 0;
}