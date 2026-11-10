# Aurora missing transforms inside the volume: literature and documentation

Reviewed2026-11-10. User asks whether literature discusses loss when the arm
moves toward the field generator. The sources below support pose/orientation-
dependent or interference-related loss, including inside the working volume.
They do not establish a universal rule that approaching the generator causes
loss, or uniquely diagnose this project's Planar2.0 setup.

## Closest direct experimental match

Koushik Kanti Mandal(2016), *Intra-Operative Needle Tracking Using Optical
Shape Sensing Technology*, master's thesis, Polytechnique Montréal.
Chapter6,page84,Figure6.2. Reports Aurora failing to track certain needle
orientations relative to the generator while the tip is inside the working
volume; the tracker screenshot shows BadFit. An approximately perpendicular
orientation is one failed example. This is evidence of in-volume tracking
failure, not a proof of the same mechanism in our sensor/generator.

-Institutional record:https://publications.polymtl.ca/2090/
-PDF:https://publications.polymtl.ca/2090/1/2016_KoushikKantiMandal.pdf
  Note:84is the printed thesis page; PDF-viewer indexing includes front matter
  and can differ. SearchFigure6.2 rather than relying on the URLpage anchor.

The PDF fetch timed out in direct webopen, but the indexed full-text find
returned the relevant passage and figure caption; its list of figures assigns
printedpage84. No local full PDF was claimed to have been downloaded.

## Continuum-robot measurement interruptions

Lachlan Scott,Tangyou Liu,Liao Wu(2024), *A Low-Cost Teleoperable Surgical Robot
with a Macro-Micro Structure and a Continuum Tip for Open-Source Research*,
IEEE/ASMEAIM2024,pp847-852,DOI10.1109/AIM55361.2024.10636985.
PDFpage5 discusses brief Aurora BadFit interruptions and sparse usable
measurement data. They attempted a~10cm tip-to-generator spacing and a
perpendicular needle-axis orientation. This is an attempted mitigation,
not a demonstrated universal cure. Different tool geometry/setup can make
orientation recommendations differ from Mandal's failed example.

-Authorpreprint:https://arxiv.org/abs/2405.16084
-PDF:https://arxiv.org/pdf/2405.16084#page=5
-ProceedingsPDF:https://ras.papercept.net/images/temp/AIM/files/0264.pdf

## NDI manufacturer guidance

*Aurora V2 User Guide*,sections5.11 and12, manufacturer-authored document
available through a third-party mirror:
https://www.scribd.com/document/1010309716/Aurora-V2-User-Guide

Missing is failure to calculate a transformation, not synonymous with an
out-of-volume event. Listed causes include excessive motion, interference,
system errors and detection-region limits. Troubleshooting also addresses
tools disappearing before nominal volume edges and environmental distortion.
Cable instructions give30mm minimum separation between tool cables and the
FieldGenerator cable, especially avoiding parallel runs. This is cable-to-
cable clearance, NOT a30mm sensor-to-generator tracking distance.

This is an older V2 guide mirrored outside NDI, not a verified current
Planar2.0/V3 installation manual. Use its general diagnostics as supporting
context and confirm installed-hardware limits with the appropriate NDI guide.

## Official brochure: accuracy is environment-dependent

NDI's2019Aurora brochure states performance depends on tool design and metal,
with specifications collected in an electromagnetically undisturbed setting.
https://www.ndieurope.com/wp-content/uploads/2019/11/8300163_rev009_Aurora.pdf
This supports checking environmental distortion but does not demonstrate
near-generator dropout in our hardware. Its planar drawing includes a50mm
offset; do not transfer that value uncritically between generator models.

## Interpretation for our experiments

Recorded tip losses occur with clear volume flags, without valve actuation,
and transfer fromP3toP2after nominal180degree arm rotation. Those are our own
observations. The literature supports treating them as tracker-relative
pose/orientation/quality failures, rather than concluding that a fixedP3valve
is responsible or the sensor must be outside the nominal volume.

Moving toward the generator changes both distance and sensor orientation in
the current arm trajectory. The experiments have not separated those two
variables. A useful diagnostic comparison is orientation changes at nearly
fixed position, then position/distance changes at nearly fixed orientation,
as the physical setup permits. Holding an intermediate failing pose is also
more discriminating than only holding the large-bend endpoint. These are
proposals, not hardware actions or validated fixes.

No universal failed angle, distance threshold, firmware defect or sensor
damage is inferred from the cited work. No tracking gates, runtime code or
hardware settings were changed. No device/model access,commit orpush.
