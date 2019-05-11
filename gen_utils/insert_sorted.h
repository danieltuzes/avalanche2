#include "stdafx.h"

template <typename T, typename Comp>
void insert_sorted(std::vector<T>& vec, T const& e, Comp const& comp) {
	auto const it = std::lower_bound(vec.begin(), vec.end(), e, comp);

	if (it != vec.end() && !comp(e, *it)) { return; }

	vec.insert(it, e);
}

template <typename T>
void insert_sorted(std::vector<T>& vec, T const& e) {
	insert_sorted(vec, e, std::less<T>{});
}