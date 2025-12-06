# AutoC / crrcsim Coordinate & Control Conventions

This note captures the shared assumptions between the AutoC GP framework and the
crrcsim‐based evaluation environment. It is intended as a quick reference when
analysing logs or integrating new flight data.

## Frames & Axes
- **World frame:** North–East–Down (NED).  
  - `+X` points toward geographic north.  
  - `+Y` points toward geographic east.  
  - `+Z` points downward toward the centre of the earth.
- **Body frame:** Right-handed aircraft axes.  
  - `+x_b` forward through the nose.  
  - `+y_b` out the right wing.  
  - `+z_b` downward through the belly.

## Position & Velocity Units
- crrcsim publishes FDM state in feet / feet-per-second.  
- The AUTOC bridge converts to metres / metres-per-second before storing in
  `AircraftState`.  
- Path generator targets and GP sensors operate strictly in metres and radians.

## Orientation Representation
- **Canonical state:** Unit quaternion `(w, x, y, z)` carried through the GP
  stack and exchanged with the simulator.  
- **Euler read-out:** When needed for logging, the quaternion is converted to
  yaw–pitch–roll via Eigen’s `eulerAngles(2, 1, 0)` (Z–Y–X sequence) and wrapped
  to ±π radians.  
- **Body attitude sense:**  
  - Positive roll = right wing down.  
  - Positive pitch = nose up.  
  - Positive yaw = nose right.

## Sensor & Derived Quantities
- `GETALPHA` and `GETBETA` rotate the velocity vector into the body frame using
  the quaternion and compute angle-of-attack / sideslip with the standard
  definitions (radians).  
- `GETROLL_RAD`, `GETPITCH_RAD` expose the wrapped Euler angles in radians.  
- Distance-related primitives (`GETDHOME`, `GETDTARGET`, etc.) operate on metre
  vectors in the NED frame.

## Control Command Polarity & Scaling
- GP program outputs are clamped in the range **−1 … 1** for pitch, roll,
  throttle.  
- Before passing to crrcsim (`TSimInputs`), the bridge rescales to simulator
  stick conventions:
  - Elevator = **−pitch / 2** ⇒ matches crrcsim’s `−0.5 … 0.5`, positive value
    = trailing edge down (nose up).  
  - Aileron  = **roll / 2** ⇒ `−0.5 … 0.5`, positive value = right wing down.  
  - Throttle = **(throttle / 2) + 0.5** ⇒ `0 … 1`, positive GP command =
    increasing power.  
- Rudder, flaps, spoilers, etc. remain neutral in the current integration.

## Summary Table

| Aspect                    | AutoC GP side                        | crrcsim / FDM side                 | Notes                                     |
|---------------------------|--------------------------------------|------------------------------------|-------------------------------------------|
| World axes                | NED (m, m/s)                         | NED (ft, ft/s)                     | bridge converts ft → m                    |
| Body axes                 | Right-handed aircraft                | Right-handed aircraft              | x forward, y right, z down                |
| Orientation               | Quaternion `(w,x,y,z)`               | Quaternion from EOM01              | logging uses wrapped yaw/pitch/roll [rad] |
| Pitch command             | `−1 … 1`, positive = pull up         | Elevator `−0.5 … 0.5`, trailing edge down | mapped via `−pitch/2`              |
| Roll command              | `−1 … 1`, positive = right wing down | Aileron `−0.5 … 0.5`               | mapped via `roll/2`                       |
| Throttle command          | `−1 … 1`, positive = more power      | `0 … 1`                            | mapped via `(throttle/2)+0.5`             |
| Angle units               | Radians                              | Radians                            | derived sensors wrap to ±π                |
| Distances / positions     | Metres                               | Feet                               | NED frame                                 |
| Velocity                  | m/s (body/world)                     | ft/s                               | conversion done before GP access          |

Keep this file updated if additional actuators, sensors, or frames are added.
