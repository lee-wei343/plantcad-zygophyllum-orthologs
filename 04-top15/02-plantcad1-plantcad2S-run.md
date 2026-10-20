# PlantCAD scoring (manual step, LLMPipeGUI)

Run outside R in LLMPipeGUI. Inputs come from `01-select-top15.R`, and
`03-plot-top15.R` reads the results.

- Tool: LLMPipeGUI from <https://github.com/SergiuM002/LLMPipeGUI>
- Window: 8192 bp each
- Server: julio.genetik.uni-koeln.de
- Date run: 2026-10-01 to 2026-10-04
- Inputs: `out/tx_a_thaliana.fa`, `out/tx_z_stapffii.fa`
- Models / output folders:
  - PlantCAD1 → `out/{a_thaliana,z_stapffii}_top15_PCAD1_results/`
  - PlantCAD2 small → `out/{a_thaliana,z_stapffii}_top15_PCAD2S_results/`
