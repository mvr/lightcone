Lightcone
=========

Lightcone is a Game of Life search tool for conduits and oscillators,
in the style of `ptbsearch`, `catalyst` and CatForce.

For some details on how it works, see this [blog
post](https://mvr.github.io/posts/lightcone.html).

Usage
-----
To build, run `git submodule update --init --recursive` and then:

```
cmake .
./Lightcone examples/minimal.toml
```

Input Parameters
----------------

See `examples/` for complete working inputs.

| Parameter                | Format                  | Default    | Meaning                                                 |
|--------------------------|-------------------------|------------|---------------------------------------------------------|
| `pattern`                | `"""LifeHistory RLE"""` | (required) | Initial pattern                                         |
| `pattern-center`         | `[x, y]`                | `[0, 0]`   | Shift applied to `pattern`                              |
| `max-catalysts`          | `n`                     | `100`      | Maximum number of catalyst                              |
| `max-transparent`        | `n`                     | `0`        | Maximum number of transparent catalysts                 |
| `max-limited`            | `n`                     | `1`        | Maximum number of 'limited' catalysts placements        |
| `first-active-range`     | `[min, max]`            | `[0, 100]` | Allowed generation range for first interaction          |
| `active-window-range`    | `[min, max]`            | `[0, 500]` | Allowed duration for overall reaction                   |
| `min-stable-time`        | `n`                     | `8`        | Required recovery duration to count as a solution       |
| `continue-after-success` | `true` / `false`        | `true`     | Keep searching after finding solution?                  |
| `max-stationary-time`    | `n`                     | `0`        | Maximum generations a non-catalyst cell may sit still   |
| `max-stationary-count`   | `n`                     | `0`        | Total population that can break the previous constraint |
| `use-bloom-filter`       | `true` / `false`        | `false`    | Enable Bloom filter                                     |
| `print-summary`          | `true` / `false`        | `false`    | Print all solutions as a row-packed summary at the end  |


### Catalysts

Add one `[[catalyst]]` block per catalyst type.

| Parameter        | Format                     | Default    | Meaning                                            |
|------------------|----------------------------|------------|----------------------------------------------------|
| `rle`            | `"RLE"`                    | (required) | Catalyst pattern                                   |
| `recovery-range` | `[min, max]`               | `[0, 100]` | Allowed interaction duration                       |
| `transparent`    | `true` / `false`           | `false`    | Catalyst is transparent?                           |
| `limited`        | `true` / `false`           | `false`    | Catalyst is limited for purposes of `max-limited`? |
| `required`       | `"LifeHistory RLE"`        | none       | Cells that must remain unchanged                   |
| `approach`       | `"LifeHistory RLE"`        | none       | Single approach signature OR                       |
| `approaches`     | `["LifeHistory RLE", ...]` | none       | Multiple approach signatures                       |
| `forbidden`      | `["LifeHistory RLE", ...]` | none       | Forbidden approach patterns at placement time      |

Notes:
- For non-transparent catalysts, you should provide `required` and at least one approach (`approach` or `approaches`).
- `forbidden` can be used to prevent glider eating reactions, which can bog down a search.
- Examples also include metadata keys like `summary` for debugging purposes; unknown keys are ignored by the parser.

### Filters

Filters are checked when reporting solutions, and don't prune the
search otherwise. There are two additional global parameters:

| Parameter                    | Format            | Default | Meaning                                                     |
|------------------------------|-------------------|---------|-------------------------------------------------------------|
| `filter-mode`                | `"ALL"` / `"ANY"` | `"ALL"` | How multiple `[[filter]]` blocks combine                    |
| `filter-match-survival-time` | `n`               | `0`     | How long a `MATCH` filter must survive undisturbed to count |

And each `[[filter]]` block has the following options:

| Parameter      | Format                           | Default    | Meaning                                   |
|----------------|----------------------------------|------------|-------------------------------------------|
| `filter`       | `"LifeHistory RLE"`              | `""`       | Target pattern / mask                     |
| `filter-pos`   | `[x, y]`                         | `[0, 0]`   | Translation applied to `filter`           |
| `filter-type`  | `"EXACT"` / `"EVER"` / `"MATCH"` | `"EVER"`   | Filter type                               |
| `filter-gen`   | `n`                              | `-1`       | Single generation to check for `EXACT` OR |
| `filter-range` | `[min, max]`                     | `[-1, -1]` | Generation range to check for `EXACT`     |

Semantics for `filter-type`:
- `EXACT`: equality check during the specified generation/range
- `EVER`: equality check at any generation up to the filter horizon
- `MATCH`: match any location/symmetry during the specified generation/range.
