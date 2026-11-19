Download the ZIP in this folder and extract it first. Its internal filenames are preserved. Open MATLAB in the extracted folder, with all six main .m files beside Elcentro.mat. The Q3 alternative is separate and optional, as explained below.

# Dynamics code review

Run the .m files from the same folder as Elcentro.mat. Keep the capital E in the filename.

Order from the supplied Pastebin:
1. Q1Assignment1.m and Q1Assignment1.pdf
2. Q2Assignment1.m and Q2Assignment1.pdf
3. Q3Assignment1.m and Q3Assignment1.pdf
4. FreeDampedSolution.m and FreeDampedSolution.pdf (z = 0.20)
5. FreeDampedSolution_z025.m and FreeDampedSolution_z025.pdf (z = 0.25)
6. FreeDampedSolution2.m and FreeDampedSolution2.pdf (z = 0.40)

Each program has its own PDF. Originals in the linked repository were not changed.

Changes to review:
- Basic scripts and loops, no cumtrapz or gradient calls. Standard load/plot/zeros/sqrt/exp/sin/cos commands remain.
- Q1 calculates true relative velocity and absolute acceleration peaks instead of the reference's pseudo-spectra. Periods are 0.05 to 4 s; smaller internal steps are used for short periods.
- Q2 uses consistent column vectors, the exact data filename, the total-displacement Newmark formula, explicit RK4 slopes and explicit trapezoidal integration. Its plots compare methods directly for each mass case.
- The second-last free-vibration program changes z to 0.25 only.
- The last program divides by wd to satisfy its stated initial velocity of 8 m/s. This corrects the reference formula.
- Q3 uses the standard harmonic-load interpretation of DAF, following the supplied reference. The photo links it to the earthquake case but does not define its normalization. Confirm whether the instructor instead wants peak earthquake displacement divided by the static displacement under peak effective earthquake force. This interpretation remains unresolved.

Data assumption:
The supplied Elcentro.mat contains p and t, but no unit/station/component metadata. The conversion p*0.01 is retained from the reference, assuming cm/s^2. The provenance and units need confirming against the original record before final submission. No new PEER download is claimed.

Checks:
All six .m files ran to completion in GNU Octave 6.4 with figure creation. Saved array shapes, finite values and equation-of-motion residuals for Q2 were checked. Free-vibration arrays agree with the stated initial-condition formula. The PDFs' figures are drawn from those actual saved results. MATLAB itself was not available for an additional run.

References supplied for review:
https://github.com/shriwedsshibu-png/dynamics-Assignment
https://pastebin.com/k9knARTn
Previous coding-style examples:
https://github.com/shriwedsshibu-png/ACM-IITD
Assignment requirements: the two notebook photos supplied in the conversation.
These are review materials. Check course rules, acknowledge assistance as required, and make sure you can explain every step. No assignment submission has been made. These separate review copies were uploaded with permission; original repository files remain unchanged.

Separate Q3 alternative:
The Q3 alternative - review only folder contains an optional earthquake-normalized definition, not a replacement for Q3Assignment1. Confirm the intended definition first. To run it, copy its .m into this parent folder beside Q2Assignment1.m and Elcentro.mat.
