# Hardware-in-the-loop residual model learning: proposed extension

2026-11-03. User asked whether the local HIL soft-robot paper could improve
model accuracy. This is a reviewed proposal only, not authorization to automate
hardware, change controllers/Simulink, replace the current comparison plan or
start a new training campaign.

## Reviewed source

Jesus Marquez, Charles Sullivan, Ryan M. Price and Robert C. Roberts,
Hardware-in-the-Loop Soft Robotic Testing Framework Using an Actor-Critic
Deep Reinforcement Learning Algorithm, IEEE RA-L8(9),6076-6082,2023.
DOI10.1109/LRA.2023.3301215.
Local source: C:\Users\dperera\OneDrive - Texas A&M University\Lab\Research\fDyn_ExpCom\hyDyn_Ctrl\Hardware-in-the-Loop_Soft_Robotic_Testing_Framework_Using_an_Actor-Critic_Deep_Reinforcement_Learning_Algorithm.pdf
SHA256: 3629dbddacd75bea0e43afac75ba37839fe4f2a0926d4734fd934d8d8e87ce22
Read full extracted paper; visually reviewed pp6078-6079 (PDF pages3-4),
especially SectionII.E, SectionIII and Fig3. Timing/repeatability checked on
p6080; exploration/replay discussion on p6081.

## What the paper demonstrates

ROS-linked physical two-actuator pneu-net locomotion with model-free TD3.
Actor outputs two0..100% actuation commands mapped to syringe-plunger positions;
it does not predict residual forces or next-state dynamics. Flex readings
provide observations and AprilTag displacement supplies task reward. Replay,
noisy exploration and episodic interaction update actor/critics. Reported
improvements are locomotion-policy outcomes versus random actuation, not
identified-model accuracy or PAC/PDFB tracking results.

Physical learning was costly:3000 steps averaged about5h per training session;
9/10 repeat sessions converged. One robot exploited the reward at the expense
of distance traveled. These specifics motivate careful objective/validation
choices; they are not convergence or timing estimates for our continuum arm.

## Adaptation to this project

Borrow physical experiment -> logged transitions -> replay -> candidate fit ->
validation -> next experiment. Keep PAC/PDFB as controllers. Train the generalized-
force residual MODEL using supervised measured next-length and free-path errors
through the existing nominal dynamics/integrator, with force magnitude/rate
regularization. The GRU still outputs two generalized forces in N on the RHS.
Do not turn it into a learned pressure policy or train solely on reference
tracking reward. Reference position is not the measured next-state label.
No direct physical-force target exists in current logs; prediction evidence
alone is not measured force accuracy.

Proposed staged loop:
1. Run an operator-approved bounded trajectory with fixed neural weights and
   stable controller settings. Log measured q/estimated dq, commanded pressure,
   measured pressure if available, timestamps, validity, references, limits and
   all controller/adaptation diagnostics.
2. Add valid chronological sequences to an old+new replay dataset. For our GRU,
   sample contiguous windows with correct causal warmup and gap/state resets;
   do not copy TD3's arbitrary independent-transition sampling into recurrence.
3. Fit a candidate off the real-time control path, initially between trials.
   Use a separate worker/process only after the manual loop is reliable.
4. Evaluate on fixed held-out trials, one-step/.5-1s paths and the known weak
   circle/loading/unloading examples; also check force variation, bounded inverse
   state purity and timing. Validation recordings must not enter replay training.
5. Apply an accepted candidate between trials with an explicit measured-prefix
   recurrent-state initialization/reset. Preserve the incumbent, source/data
   hashes and candidate identity. Verify physical performance in fresh repeats.
6. Choose subsequent bounded experiments to cover observed state/pressure/speed
   weaknesses. Prediction error can guide selection; the current single GRU
   does not supply calibrated uncertainty. Learned active exploration/ensembles
   would be additional future work, not a present capability.

Use the existing MATLAB/Simulink logging and Python/MATLAB training interface as
a starting point; a ROS migration and actor-critic policy are not prerequisites.
The first controller comparison should still freeze model weights, match gains
and pressure allocation, and compare model-off/on within PAC and PDFB. Resolve
PAC adaptive compensation overlap and nominal-model consistency first.
Current GRU128 is a short-horizon candidate; GRU64 remains the reference because
long circle accuracy regressed with128. HIL data collection might improve
coverage of actual controller operating conditions, but cannot guarantee better
model accuracy or tracking. Preserve a fresh evaluation cohort before updating.
