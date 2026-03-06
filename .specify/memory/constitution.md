<!--
  ============================================================================
  SYNC IMPACT REPORT
  ============================================================================
  Version change: 1.0.0 → 1.1.0

  Modified principles: None

  Added sections:
    - Multi-repo build system documentation (GP, CRRCSim, xiao-gp)

  Removed sections: None

  Templates requiring updates:
    - .specify/templates/plan-template.md: ✅ Compatible (has Constitution Check section)
    - .specify/templates/spec-template.md: ✅ Compatible (no constitution-specific refs)
    - .specify/templates/tasks-template.md: ✅ Compatible (supports testing workflow)

  Follow-up TODOs: None
  ============================================================================
-->

# GPC++ / AutoC Constitution

## Core Principles

### I. Testing-First

All significant changes MUST include tests that validate the intended behavior. The test-driven
development workflow is:

1. Write tests that specify expected behavior before implementation
2. Verify tests fail (confirming they test something meaningful)
3. Implement the feature or fix
4. Verify tests pass
5. Refactor if needed while maintaining passing tests

**Rationale**: The GP/AutoC system involves complex mathematical operations and control logic
where subtle bugs can produce invalid evolution results. Tests provide confidence that changes
do not regress existing functionality.

**Exemptions**: Exploratory prototyping and research spikes may skip tests, but any code
promoted to mainline MUST have corresponding tests.

### II. Build Stability

All commits to the main branch MUST compile successfully and pass the test suite. The build
MUST NOT be left in a broken state across any affected repository.

**Build verification requirements**:

- GP: `cd ~/GP && make` MUST compile without errors
- CRRCSim: `cd ~/crsim/crrcsim-0.9.13/build && make` MUST compile without errors
- xiao-gp: `cd ~/xiao-gp && pio run` MUST compile without errors
- All existing tests MUST pass

**Rationale**: A broken build blocks all other development work. The multi-repo system
requires all build paths to remain functional.

**Recovery**: If a commit breaks the build, the fix is highest priority. No other work
proceeds until build stability is restored.

### III. Dual-Mode Parity

The GP tree evaluation mode and bytecode interpreter mode MUST produce identical results
for the same input conditions. Any change to evaluation logic MUST be reflected in both:

- `autoc-eval.cc` (MyGene::evaluate for GP tree traversal)
- `gp_bytecode.cc` (GPBytecodeInterpreter::evaluate for bytecode interpretation)

**Verification requirements**:

- New operators MUST be implemented in both evaluation modes
- Changes to existing operators MUST be applied to both modes
- EvaluateMode=1 verification MUST confirm identical fitness values

**Rationale**: The bytecode interpreter enables deployment of evolved controllers independent
of the full GP system. Parity ensures that deployed behavior matches evolved behavior.

## Build & Deployment

This is a **multi-repo build system**. Each repository has its own build approach.
These build instructions MUST always be followed.

### Repository Overview

| Repository | Build System | Location |
|------------|--------------|----------|
| GP (core + autoc) | CMake + Make | `~/GP` |
| CRRCSim | CMake + Make | `~/crsim/crrcsim-0.9.13` |
| xiao-gp | PlatformIO | `~/xiao-gp` |

### From-Scratch / CMakeLists.txt Changes

When building from scratch OR when `CMakeLists.txt` changes, use the rebuild scripts:

**GP/AutoC**:
```bash
# Debug mode (better debugging, faster compile)
cd ~/GP/autoc && bash rebuild.sh

# Perf mode (optimized -O3, for long-running large tests)
cd ~/GP/autoc && bash rebuild-perf.sh
```

**CRRCSim**:
```bash
cd ~/crsim/crrcsim-0.9.13
mkdir -p build && cd build
cmake ..
make
```

### Incremental Builds

For day-to-day development after initial setup, use incremental builds:

```bash
# GP (includes autoc cmake subdir)
cd ~/GP && make

# CRRCSim
cd ~/crsim/crrcsim-0.9.13/build && make

# xiao-gp (with platform knobs as needed)
cd ~/xiao-gp && pio run
```

### Build Mode Selection

| Mode | Script | Use Case |
|------|--------|----------|
| Debug | `rebuild.sh` | Development, debugging, fast iteration |
| Perf | `rebuild-perf.sh` | Long-running evolution, large population tests |

### Deployment Modes

1. **Evolution Mode** (EvaluateMode=0): Full GP tree evaluation during population evolution
2. **Verification Mode** (EvaluateMode=1): Bytecode interpreter with parity checking
3. **Standalone Bytecode**: gpextractor output for embedded deployment (xiao-gp)

### Release Checklist

- [ ] All tests pass
- [ ] GP incremental build succeeds (`cd ~/GP && make`)
- [ ] CRRCSim build succeeds (`cd ~/crsim/crrcsim-0.9.13/build && make`)
- [ ] xiao-gp build succeeds (`cd ~/xiao-gp && pio run`)
- [ ] Dual-mode parity verified (EvaluateMode=0 and EvaluateMode=1 produce same results)
- [ ] CLAUDE.md updated if build/config changes

## Governance

### Amendment Process

1. Proposed amendments MUST be documented with rationale
2. Amendments affecting core principles require explicit approval
3. All dependent artifacts (specs, plans, tasks) MUST be reviewed for consistency

### Versioning Policy

The constitution follows semantic versioning:

- **MAJOR**: Removal or redefinition of core principles
- **MINOR**: New principle or section added; material guidance expansion
- **PATCH**: Wording clarifications, typo fixes, non-semantic refinements

### Compliance Review

- All PRs MUST verify compliance with applicable principles
- Constitution violations MUST be documented and justified if exempted
- The plan-template.md "Constitution Check" section enforces review

### Runtime Guidance

For day-to-day development guidance, refer to [CLAUDE.md](../../CLAUDE.md) which contains
build commands, architecture details, and common workflows.

**Version**: 1.1.0 | **Ratified**: 2026-03-06 | **Last Amended**: 2026-03-06
