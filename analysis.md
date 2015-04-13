1. Pick the **JH** gene with the highest score and a minimum of 30bp.

2. Choice of **JH** sets the reading frame for the sample.

3. Pick the **VH** gene with the highest score and a minimum of 40bp and 70% identity. Do we want to filter on *functionality* or *in-frame*?

4. Penalize **DH** scores for overlaps with **VH** and **JH** and trim the hits of the overlap. Penalty is 0.5 * **overlap**.

5. Pick a **DH** gene with the highest score that has at least 4bp match in a row.

6. Check if the chosen **VH** and **JH** are in frame

7. Check for stop codons encoded in sample.

8. Annotate *hybrid* samples. Does than mean more than 1, non overlapping, **DH** matches?

9. Identify CDR3. Scan V region backwards for Cycstine. Scan J region for **TGG GG** in frame. 

10. Annotate palindromic sequences on the ends of genes.
