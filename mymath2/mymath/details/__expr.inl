#pragma once

namespace mymath {
	namespace {
		template<class T, class U>
		constexpr const T& min(const T& a, const U& b) {
			return a < b ? a : b;
		}

		template<class T, class U>
		constexpr const T& max(const T& a, const U& b) noexcept {
			return a > b ? a : b;
		}
		 
	}
}