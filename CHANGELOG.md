# Changelog


## [0.2.0] - 2018-10-26

### Changed
- Tasks are generated in place on each rank.
- UVEngines are reused by the worker loops.
- Alt/Az quantities are only recalculated for each source when the time changes.
- Use `reuse_splines` option on UVBeams to reduce interpolation time.
- Beams are scattered as strings, not objects.


## [0.1.0] - 2018-10-24
