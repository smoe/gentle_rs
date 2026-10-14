# Recorded process and timing observations

These values are transcribed from the original recorded tool calls and messages.
No GUI process was restarted for this supplement.

## Microarray

- Input: retained public 83,686-bp locus project; disposable copy hash is in
  `../native-x11/microarray/input.sha256`.
- Main project window opened and `01-overview.png` was captured in the original
  evidence (the image is deliberately omitted here).
- After opening the sequence window, a 3-second `xwininfo -root` probe timed out
  with exit code 124.
- First process sample: GENtle 48.1% CPU, 1,281,072 KiB RSS, state `Sl+`.
- At elapsed 131 seconds: GENtle 36.0% CPU, 1,281,072 KiB RSS, state `Sl+`;
  a 5-second X11 probe again timed out with exit code 124.
- The bounded attempt was stopped. No source state was overwritten. Native
  navigation/results/save-reopen and readable result capture remain open.

## Gel

- Input: retained synthetic project with two embedded gel images; disposable
  copy hash is in `../native-x11/gel/input.sha256`.
- The project loaded and the genuine `Patterns -> Gel Image Analysis...` route
  was exercised.
- At elapsed 114 seconds: GENtle 34.1% CPU, 412,176 KiB RSS, state `Sl+`.
- A 3-second `xwininfo -root` probe timed out with exit code 124.
- The editor still had not become X11-responsive after more than two minutes;
  the bounded attempt was stopped. Calibrated/null result display and GUI
  save/reopen remain open.

These are environment-bounded presentation observations, not generalized
product-defect claims.
