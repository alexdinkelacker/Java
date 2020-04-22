package dnasplicing;

public class DnaNode {
	public String dnaSequence;
	public DnaNode previous;
	public DnaNode next;


	public DnaNode(String initialDnaSequence) {
		dnaSequence = initialDnaSequence;
	}
}
