.. This document is written in reStructuredText.
.. Build command:
   rst2html.py --date --stylesheet-path ./style.css readme.rst readme.html

======================================================
PAPA - a program for predicting prion-forming proteins
======================================================

--------------------------------------------------------

PAPA predicts the prion-forming propensity of a protein using a
sliding-window approach.
Using a set of known prion-forming and non-prion-forming domains the
prion propensity of each amino acid was established.
Each position in a protein is characterized by the average prion
propensity in a length fourty-one window around it.
Prediction scores are computed by averaging the scores of fourty-one
consecutive windows, and the maximum across the protein is its
prion-forming propensity.
Segments with a score above 0.05 have a high probability of forming
prion-like aggregates.  The default setting only scores segments that
are predicted to be disordered according to `FoldIndex <http://bip.weizmann.ac.il/fldbin/findex>`_; if no such segment
is found, PAPA provides a prion propensity score of -1.  This disorder
requirement can be turned off.  
Please note that PAPA has only been validated for Q/N-rich sequences.



Download
------------

The papa program is a command-line Python script available 
`here <./papa.tgz>`_.
To test the program you can download 
`a few of the sequences <./sequences.fasta>`_
mentioned in the paper.

Usage
-----

Typical papa usage:

   ``python papa.py -o results_file fasta_file``

where:

fasta_file: 
    path to a fasta-format file containing the protein sequences
results_file:
    path to a file containing the output
    The output is a comma delimited file whose columns are:
    sequence id, maximum score, position.
    Maximum score is the largest window score, and position is the position where it occurs.

For full usage type: ``python papa.py -h``.  Note that windows users
may need to provide the full 
   
License
-------

GPLv3_

.. _GPLv3: http://gplv3.fsf.org/

All programs in this collection are free software: 
you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Acknowledgements
----------------

This material is based upon work supported by the National Science Foundation under Grant No. 1023771.  Any opinions, findings and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation (NSF).

