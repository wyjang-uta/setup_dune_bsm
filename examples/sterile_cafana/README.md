<h1>DUNE BSM sterile analysis (with CAFAna) examples</h1>

<p>
NOTE: The required framework tools are not yet in master; need to use feature/wallbank_SterileAna branch.

Two macros are currently included, which can be used to make and analyze oscillated spectra.  Run as: <br>
cafe -bq -l <n> MakeOscSpectra.C <br>
cafe -bq AnalyzeOscSpectra.C <file_name>.root <br>
where n is the number of input files to run over, and file_name is the output file produced by the first macro.
</p>