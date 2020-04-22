package sortcomparison;

import java.util.*;

public class BasicSorter implements Sorter {

	/**
	 * Uses the insertion sort algorithm to sort part or all of a list in ascending
	 * order.
	 * 
	 * @param data
	 *            The list of elements to sort
	 * @param fi
	 *            Index of the first element to be sorted.
	 * @param n
	 *            The number of elements in the section to be sorted.
	 */
	@Override
	public void insertionSort(String[] data, int fi, int n) {
		for (int i = fi; i < fi + n; i++) {
			for (int j = i; j > fi && data[j].compareTo(data[j - 1]) < 0; j--) {
				swap(data, j, j - 1);
			}
		}
	}


	private void quickInsertionSort(String[] data, int fi, int n) {
		for (int i = fi + 1; i < n; i++) {
			String key = data[i];
			int j = i - 1;

			while (j >= fi && data[j].compareTo(key) > 0) {
				// data[j + 1] = data[j];
				System.arraycopy(data, j, data, j + 1, 1);
				j--;
			}
			data[j + 1] = key;
		}
	}


	/**
	 * Uses the quick sort algorithm to sort all or part of the data in ascending
	 * order. quickSort() must:
	 * 
	 * a) Call the insertionSort() function for segments of size 15 or less.
	 * 
	 * b) Call the partition() function to do its partitioning.
	 * 
	 * @param data
	 *            The list of elements to sort
	 * @param fi
	 *            Index of the first element to be sorted.
	 * @param n
	 *            The number of elements in the section to be sorted.
	 */
	@Override
	public void quickSort(String[] data, int fi, int n) {

		quickSorter(data, fi, n - 1);

	}


	private void quickSorter(String[] data, int fi, int n) {
		if (fi >= n) {
			return;
		}
		if (n - fi <= 15) {
			quickInsertionSort(data, fi, n);
		}

		// partition
		int i = fi, j = n;
		String pivot = getMedian(data, fi, n);

		while (i <= j) {
			while (data[i].compareTo(pivot) < 0) {
				++i;
			}
			while (data[j].compareTo(pivot) > 0) {
				--j;
			}

			if (i <= j) {
				String tmp = data[i];
				data[i++] = data[j];
				data[j--] = tmp;
			}
		}

		quickSorter(data, i, n);
		quickSorter(data, fi, j);
	}


	private String getMedian(String[] data, int left, int right) {
		int center = left + (right - left) / 2;

		if (data[left].compareTo(data[center]) > 0) {
			swap(data, left, center);
		}

		if (data[left].compareTo(data[right]) > 0) {
			swap(data, left, right);
		}

		if (data[center].compareTo(data[right]) > 0) {
			swap(data, center, right);
		}

		swap(data, center, right);
		return data[right];
	}


	/**
	 * Partitions part (or all) of the list. The range of indices included in the
	 * partitioning is [firstIndex, firstIndex + numberToPartition -1].
	 * 
	 * partition() ** must ** use the median-of-three method to prevent O(n^2)
	 * running time for sorted data sets.
	 * 
	 * 
	 * 
	 * @param data
	 * @param fi
	 * @param n
	 * @return The index, relative to data[0], where the pivot value is located at
	 *         the end of this partitioning.
	 */
	@Override
	public int partition(String[] data, int fi, int n) {
		int left = fi;
		int right = fi + n - 1;
		int middle = (fi + (n - 1)) / 2;

		if (data[left].compareTo(data[middle]) > 0) {
			swap(data, left, middle);
		}
		if (data[left].compareTo(data[right]) > 0) {
			swap(data, left, right);

		}
		if (data[middle].compareTo(data[right]) > 0) {
			swap(data, middle, right);
		}

		swap(data, middle, right);

		String pivot = data[right];

		int i = left;
		for (int j = left; j < right; j++) {
			if (data[j].compareTo(pivot) <= 0) {
				swap(data, i, j);
				i++;
			}
		}

		swap(data, i, right);
		return i;

	}


	private void swap(String[] data, int i, int j) {

		String temp = data[i];
		data[i] = data[j];
		data[j] = temp;

	}


	/**
	 * Uses the merge sort algorithm to sorts part or all of a list in ascending
	 * order. mergeSort() must:
	 * 
	 * a) Call the insertionSort() segments of size 15 or less.
	 * 
	 * b) Call the merge() function to do its merging.
	 * 
	 * @param data
	 *            list of elements to sort
	 * @param fi
	 *            Index of the first element to be sorted.
	 * @param n
	 *            The number of elements in the section to be sorted.
	 */
	@Override
	public void mergeSort(String[] data, int fi, int n) {

		int end = fi + n - 1;

		if (end <= fi || end >= data.length) {
			return;
		}

		if (n > 15) {
			int middle = (fi + end) / 2;
			mergeSort(data, fi, middle - fi + 1);
			mergeSort(data, middle + 1, end - middle);
			merge(data, fi, middle - fi + 1, end - middle);
		} else {
			insertionSort(data, fi, n);
		}

	}


	/**
	 * Merges two sorted segments into a single sorted segment. The left and right
	 * segments are contiguous.
	 * 
	 * @param data
	 *            The list of elements to merge
	 * @param fi
	 *            Index of the first element of the left segment.
	 * @param nl
	 *            The number of elements in the left segment.
	 * @param nr
	 *            The number of elements in the right segment.
	 */
	@Override
	public void merge(String[] data, int fi, int nl, int nr) {
		String[] left = new String[nl];
		String[] right = new String[nr];

		System.arraycopy(data, fi, left, 0, nl);
		System.arraycopy(data, fi + nl, right, 0, nr);

		int i = 0;
		int j = 0;
		int k = fi;
		while (i < nl && j < nr) {
			if (left[i].compareTo(right[j]) <= 0) {
				data[k] = left[i];
				i++;
			} else {
				data[k] = right[j];
				j++;
			}
			k++;
		}

		while (i < nl) {
			// data[k] = left[i];
			System.arraycopy(left, i, data, k, 1);
			k++;
			i++;

		}
		while (j < nr) {
			// data[k] = right[j];
			System.arraycopy(right, j, data, k, 1);
			k++;
			j++;

		}
	}


	/**
	 * EXTRA CREDIT
	 * 
	 * Use the heap sort algorithm to sort the list in ascending order. This method
	 * must use the heapify() method for initial heapification.
	 * 
	 * @param data
	 *            The list of elements to sort
	 */
	@Override
	public void heapSort(String[] data) {
		// TODO Auto-generated method stub

	}


	/**
	 * EXTRA CREDIT
	 * 
	 * Heapifies the given list.
	 * 
	 * @param data
	 *            The list to heapify.
	 */
	@Override
	public void heapify(String[] data) {
		// TODO Auto-generated method stub

	}

}
