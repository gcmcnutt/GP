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
- Xiao/MSP mapping: GP −1…1 is sent as RC 1000…2000 µs. Roll is direct
  (`1500 + roll*500`), pitch is inverted (`1500 − pitch*500`), throttle is
  direct (`1500 + throttle*500`).
- Rudder, flaps, spoilers, etc. remain neutral in the current integration.

## Summary Table

| Aspect                    | AutoC GP side                        | crrcsim / FDM side                 | Notes                                     |
|---------------------------|--------------------------------------|------------------------------------|-------------------------------------------|
| World axes                | NED (m, m/s)                         | NED (ft, ft/s)                     | bridge converts ft → m                    |
| Body axes                 | Right-handed aircraft                | Right-handed aircraft              | x forward, y right, z down                |
| Orientation               | Quaternion `(w,x,y,z)`               | Quaternion from EOM01              | logging uses wrapped yaw/pitch/roll [rad] |
| Pitch command             | `−1 … 1`, positive = pull up         | Elevator `−0.5 … 0.5`               | mapped via `−pitch/2` (sim) / `1500−pitch*500` (MSP) |
| Roll command              | `−1 … 1`, positive = right wing down | Aileron `−0.5 … 0.5`                | mapped via `roll/2` (sim) / `1500+roll*500` (MSP)    |
| Throttle command          | `−1 … 1`, positive = more power      | `0 … 1` (sim) / `1000…2000 µs` (MSP) | mapped via `(thr/2)+0.5` (sim) / `1500+thr*500` (MSP) |
| Angle units               | Radians                              | Radians                            | derived sensors wrap to ±π                |
| Distances / positions     | Metres                               | Feet                               | NED frame                                 |
| Velocity                  | m/s (body/world)                     | ft/s                               | conversion done before GP access          |

Keep this file updated if additional actuators, sensors, or frames are added.

## Quaternion & Euler sign conventions (VERIFIED via bench testing 2025-12-20)

### Standard Aerospace Convention (NED/FRD)
- **World frame:** NED (North-East-Down), right-handed
- **Body frame:** FRD (Forward-Right-Down), right-handed
- **Quaternion format:** `(w, x, y, z)` representing **earth→body** rotation
- **Right-hand rule rotations:**
  - **Roll** (φ, about +X body): Positive = right wing down
  - **Pitch** (θ, about +Y body): Positive = nose up
  - **Yaw** (ψ, about +Z body): Positive = nose right (clockwise from above)

### Actual Implementation (Bench-Tested 2025-12-20)

#### INAV (xiao-gp/msplink.cpp)
- **MSP2_INAV_LOCAL_STATE** sends quaternion as `(q0, q1, q2, q3)` = `(w, x, y, z)`
- **Convention:** Body→earth in NED/FRD (CONFIRMED via bench testing)
  - Despite INAV source code suggesting earth→body, empirical tests prove body→earth
  - The `imuTransformVectorBodyToEarth` function name is misleading regarding quaternion direction
- **Transformation applied:** Full conjugate to convert body→earth to earth→body
  - `neuQuaternionToNed()` applies: `gp_quat(q[0], -q[1], -q[2], -q[3])`
  - This conjugate operation inverts the rotation direction as required by the GP contract
- **Board alignment:** INAV applies sensor→body rotation BEFORE quaternion formation
  - Current config: `align_board_roll=1700, align_board_yaw=900`
  - Compensates for IMU mounted upside-down, rotated 90° (sensor X-axis → right wing)
  - After alignment, quaternion is in aircraft body frame (not sensor frame)
- **INAV Configurator UI:** Shows **negative pitch** for nose up (cosmetic display choice only)

#### GP/CRRCSim Contract (Standard NED/FRD)
- **Required quaternion:** Earth→body in NED/FRD with standard right-hand rule
- **Validation:** All bench tests confirm proper aerospace conventions after conjugate
- **Frame consistency:** Both INAV (after transform) and CRRCSim use identical conventions

### Bench Test Results (Verification Table)

Physical craft orientations tested with INAV board alignment `1700/0/900` and conjugate transformation:

| Orientation | INAV Sends (body→earth) | After Conjugate (earth→body) | Euler (φ,θ,ψ) | Validates |
|-------------|------------------------|----------------------------|---------------|-----------|
| Level, nose north, canopy up | `[1, 0, 0, 0]` | `[1, 0, 0, 0]` | (0°, 0°, 0°) | ✅ Identity baseline |
| Nose up 45° | `[.93, 0, -.38, 0]` | `[.93, 0, .38, 0]` | (0°, +45°, 0°) | ✅ +Y = nose up (RHR) |
| Right wing down 45° | `[.92, .35, 0, 0]` | `[.92, -.35, 0, 0]` | (+45°, 0°, 0°) | ✅ +X = right wing down (RHR) |
| Nose east (yaw 90°) | `[.7, 0, 0, -.7]` | `[.7, 0, 0, .7]` | (0°, 0°, +90°) | ✅ +Z = nose right (RHR) |
| Inverted | `[0, 1, 0, 0]` | `[0, -1, 0, 0]` | (±180°, 0°, 0°) | ✅ 180° roll |
| Nose down 30°, right 20° | `[.94, -.16, .3, -.14]` | `[.94, .16, -.3, .14]` | (+20°, -30°, 0°) | ✅ Combined axes |

**Conclusion:** After conjugate transformation, quaternions follow proper NED/FRD right-hand rule and match standard aerospace conventions. The GP receives earth→body quaternions as required by the training data contract.

### Euler Angle Extraction
- Use Eigen's `eulerAngles(2, 1, 0)` for ZYX intrinsic sequence (yaw-pitch-roll)
- Wrap to ±π radians for logging
- Order: `[yaw, pitch, roll]` = `[ψ, θ, φ]`

## CRRCSim EOM01 Quaternion Convention (VERIFIED 2025-12-20)

### Mathematical Analysis

The EOM01 flight dynamics model in CRRCSim uses standard earth→body quaternion convention, confirmed by analyzing the implementation:

#### Initialization from Euler Angles (eom01.cpp:111-118)
```cpp
e_0 = cos(Psi*0.5)*cos(Theta*0.5)*cos(Phi*0.5) + sin(Psi*0.5)*sin(Theta*0.5)*sin(Phi*0.5);
e_1 = cos(Psi*0.5)*cos(Theta*0.5)*sin(Phi*0.5) - sin(Psi*0.5)*sin(Theta*0.5)*cos(Phi*0.5);
e_2 = cos(Psi*0.5)*sin(Theta*0.5)*cos(Phi*0.5) + sin(Psi*0.5)*cos(Theta*0.5)*sin(Phi*0.5);
e_3 = -cos(Psi*0.5)*sin(Theta*0.5)*sin(Phi*0.5) + sin(Psi*0.5)*cos(Theta*0.5)*cos(Phi*0.5);
```

Where: `Psi` = yaw (ψ), `Theta` = pitch (θ), `Phi` = roll (φ)

This follows the **standard ZYX (3-2-1) Euler sequence** for earth→body rotation in aerospace:
- First rotate about earth Z-axis (yaw)
- Then rotate about intermediate Y-axis (pitch)
- Finally rotate about body X-axis (roll)

The formula matches the standard aerospace convention for quaternion from Euler angles representing **earth→body** transformation.

#### Quaternion Integration (eom01.cpp:240-243)
```cpp
e_dot_0 = 0.5*( -v_R_omega_total.r[0]*e_1 - v_R_omega_total.r[1]*e_2 - v_R_omega_total.r[2]*e_3 );
e_dot_1 = 0.5*(  v_R_omega_total.r[0]*e_0 - v_R_omega_total.r[1]*e_3 + v_R_omega_total.r[2]*e_2 );
e_dot_2 = 0.5*(  v_R_omega_total.r[0]*e_3 + v_R_omega_total.r[1]*e_0 - v_R_omega_total.r[2]*e_1 );
e_dot_3 = 0.5*( -v_R_omega_total.r[0]*e_2 + v_R_omega_total.r[1]*e_1 + v_R_omega_total.r[2]*e_0 );
```

Where `v_R_omega_total` is the angular velocity in body frame (p, q, r).

These equations match the **standard quaternion kinematics for earth→body quaternion** with body-frame angular rates:
```
q̇ = 0.5 * q ⊗ ω_body
```

This is the correct formulation for integrating an earth→body quaternion using body-frame angular velocities.

#### Direct Access (eom01.h:86-89)
```cpp
virtual double getQuatW() { return e_0; }
virtual double getQuatX() { return e_1; }
virtual double getQuatY() { return e_2; }
virtual double getQuatZ() { return e_3; }
```

The quaternion is passed **directly** to the GP system with no transformations:
- EOM01 stores `(e_0, e_1, e_2, e_3)` as earth→body quaternion
- CRRCSim bridge reads via `getQuat*()` accessors (inputdev_autoc.cpp:455)
- GP receives identical quaternion in `AircraftState`

### GP Sensor Verification (gp_evaluator_portable.cc:226-240)

#### GETALPHA (Angle of Attack)
```cpp
case GETALPHA: {
    gp_vec3 velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
    result = fastAtan2(-velocity_body.z(), velocity_body.x());
    break;
}
```

**Analysis:**
- `getOrientation()` returns the earth→body quaternion
- `.inverse()` creates body→earth quaternion
- Multiplying world-frame velocity by body→earth quaternion transforms it to body frame
- This is **correct usage** of an earth→body quaternion

#### GETBETA (Sideslip Angle)
```cpp
case GETBETA: {
    gp_vec3 velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
    result = fastAtan2(velocity_body.y(), velocity_body.x());
    break;
}
```

**Analysis:**
- Uses same earth→body quaternion inversion to transform velocity to body frame
- Computes sideslip angle from body-frame velocity components
- This is **correct usage** of an earth→body quaternion

### Contract Compliance Summary

✅ **EOM01 (CRRCSim FDM):** Produces earth→body quaternions via standard ZYX Euler initialization and standard quaternion kinematics with body-frame rates

✅ **CRRCSim Bridge:** Passes quaternions directly from EOM01 to GP with no transformations

✅ **GP Sensors:** Use `.inverse()` correctly to transform world-frame vectors to body frame, confirming they expect earth→body quaternions

✅ **INAV (after conjugate):** Now produces earth→body quaternions matching the CRRCSim/GP contract

✅ **Renderer (after fix):** Correctly converts INAV blackbox body→earth quaternions to earth→body for visualization

**Conclusion:** The entire system follows a consistent earth→body quaternion convention in NED/FRD coordinates with standard aerospace right-hand rule. The INAV conjugate transformation (implemented 2025-12-20) ensures the real hardware matches the training data contract.

## Renderer Quaternion Handling (FIXED 2025-12-20)

The renderer (renderer.cc) expects earth→body quaternions in `AircraftState`, matching the GP contract:

### Visualization Usage
- **Line 1054:** `s.getOrientation() * -vec3::UnitZ()` - Rotates body-up vector to earth frame for ribbon orientation
- **Lines 2315-2319:** Attitude indicator rotates body frame axes to earth frame to extract pitch/roll

Both operations correctly use earth→body quaternions by multiplying body-frame vectors.

### INAV Blackbox Loading (renderer.cc:1453-1461)
```cpp
// INAV blackbox logs raw body->earth quaternion in NED frame
quat inavQuat(qw, qx, qy, qz);
inavQuat.normalize();

// Convert body->earth to earth->body to match renderer/GP contract
quat earthToBody(inavQuat.w(), -inavQuat.x(), -inavQuat.y(), -inavQuat.z());
earthToBody.normalize();
```

**Previous Bug:** Applied incorrect NEU→NED conversion plus conjugate, resulting in wrong orientation display.

**Fix:** Simple conjugate transformation, identical to msplink.cpp, converts INAV's body→earth to earth→body.

### CRRCSim/GP Evaluation Data Loading
CRRCSim `EvalResults` already contain earth→body quaternions (verified above), so no transformation is needed when loading evaluation data - the renderer uses them directly.
