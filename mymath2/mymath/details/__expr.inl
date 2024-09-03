#pragma once

namespace mymath {
	namespace {
		template<class T>
		constexpr const T& min(const T& a, const T& b) {
			return a < b ? a : b;
		}

		  
	}
}