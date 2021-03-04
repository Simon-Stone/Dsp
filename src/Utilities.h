#pragma once
#include <utility>
#include <vector>

namespace dsp
{
	namespace window
	{
		/// @brief Utility functions (probably only interesting for DSP-internal use)
		namespace utilities
		{
			/// @brief Extend window length by 1 sample if needed for DFT-even symmetry		
			std::pair<unsigned, bool> extend(unsigned N, bool sym);

			/// @brief Truncate window by 1 sample if needed for DFT - even symmetry
			template<class T>
			std::vector<T> truncate(std::vector<T>& w, bool needed)
			{
				if (needed)
				{
					w.pop_back();
					return w;
				}
				return w;					
			}			
		}
	}	
}
