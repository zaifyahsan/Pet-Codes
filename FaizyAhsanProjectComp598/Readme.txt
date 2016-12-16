MSA RNA Data files:
align1.txt: RF01065, 
align2.txt: RF01624, 
align3.txt: RF01068, 
align4.txt: RF00163

Sample Stockholm file: RF01065.stockholm.txt

Tree files:
align1.tree.txt, align2.tree.txt, align3.txt, align4.txt

Python file:
rnaAncestor.py

For Help:

python rnaAncestor.py -h

For Run:

Objective 1
python rnaAncestor.py RF01065.stockholm.txt align1.tree.txt False

Objective 2-8

python rnaAncestor.py align*.txt align*.tree.txt True

where * can be from 1 to 4