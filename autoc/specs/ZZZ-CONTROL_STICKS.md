# add flight controls to the renderer.cc playback system
Today we render a wallclock time in the upper right of the display.  Extend this to also play back the pitch/roll stick and the throttle control.

## Instructions
- source is in ~/GP/autoc
- examine the renderer.cc file for the playback system
- identify where in the playback records the pitch/roll/throttle are recorded
- extend the wallclock display to also show the pitch/roll stick and throttle control
- update the display in real-time as the controls change
- ensure the display is clear and readable during playback, same font and size as the wallclock
- the control stick should look like on a RC transmitter.  The control stick is x/y with center being neutral controls, and a stick left rolls left and a stick back pitches up
- the throttle control should be a vertical bar, with the bottom being zero throttle and the top being full throttle
- place the controls to the left of the stopwatch (left-to-right: stick, throttle, attitude, velocity, stopwatch) so nothing overlaps the clock; velocity/attitude must not sit on top of the clock
- when multiple arenas play back, show which arena feeds the stick/throttle HUD; default to arena 1 and switch when the user changes focus with `f` or the arrow keys
- add an attitude indicator between the throttle and the clock: blue sky on top, red ground on bottom, derived from the quaternion (no Euler shortcuts), drawn like a little horizon ball
- add a velocity indicator (m/s) between the throttle and the clock; forward speed proxy is fine (x-direction or tape/number) as long as it updates with playback
