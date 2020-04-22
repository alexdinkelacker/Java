package dnasplicing;

import static org.apache.commons.lang3.StringUtils.*;

public class LinkedDnaStrand implements DnaStrand {

	private DnaNode head = null;
	private DnaNode tail = null;
	private int nucleotideCount = 0;
	private int nodeCount = 0;
	private int appendCount = 0;


	public LinkedDnaStrand(String dnaSequence) {
		head = new DnaNode(dnaSequence);
		nucleotideCount += dnaSequence.length();
		nodeCount++;
	}


	@Override
	public long getNucleotideCount() {

		return nucleotideCount;
	}


	@Override
	public void append(String dnaSequence) {
		var newNode = new DnaNode(dnaSequence);

		DnaNode cursor = head;

		// make cursor equal the last node in the list
		while (cursor.next != null) {
			nucleotideCount += cursor.dnaSequence.length();
			cursor = cursor.next;
		}

		newNode.previous = cursor;
		cursor.next = newNode;
		newNode.next = null;
		tail = newNode;

		nucleotideCount += dnaSequence.length();
		nodeCount++;
		appendCount++;

	}


	@Override
	public DnaStrand cutSplice(String enzyme, String splicee) {
		DnaNode cursor = head;
		String sequence = head.dnaSequence;

		// First thing is scan for the enzyme
		boolean hasEnzyme = sequence.contains(enzyme);
		if (hasEnzyme) // iff this happens -- splice
		{
			boolean hasEnzymeAtTail = sequence.substring(sequence.length() - enzyme.length()).equals(enzyme);

			String[] vals = sequence.split(enzyme);
			for (String looper : vals) {
				if (!looper.isEmpty()) {
					this.append(looper);
					cursor = cursor.next;
				}
				this.append(splicee);

				nucleotideCount -= enzyme.length();
				cursor = cursor.next; // otherwise it stays one behind
			}

			if (!hasEnzymeAtTail) {
				cursor = cursor.previous;
				cursor.next = null;

				// Decrease counters for last node and first node that was removed
				appendCount--;
				nodeCount--;
				nucleotideCount -= splicee.length();
			}

			appendCount--;
			nodeCount--;
			nucleotideCount -= sequence.length();
			head = head.next;// move one forward, let the garbage collector
			// do the magic for the old head

		}

		return this;
	}


	@Override
	public DnaStrand createReversedDnaStrand() {

		LinkedDnaStrand newList = new LinkedDnaStrand(reverse(tail.dnaSequence));

		DnaNode tailCur = tail;

		while (tailCur.previous != null) {

			newList.append(reverse(tailCur.previous.dnaSequence));

			tailCur = tailCur.previous;
		}
		return newList;
	}


	@Override
	public int getAppendCount() {
		return appendCount;
	}


	@Override
	public DnaNode getFirstDnaNode() {
		return head;
	}


	@Override
	public int getDnaNodeCount() {

		return nodeCount;
	}


	@Override
	public String toString() {
		var sb = new StringBuilder();

		DnaNode cursor = head;
		while (cursor != null) {
			sb.append(cursor.dnaSequence);
			cursor = cursor.next;

		}
		return sb.toString();

	}

}
