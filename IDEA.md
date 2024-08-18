Lightcone
=========

Here's an idea for a version of CatForce that isn't so dumb. The
overhead might outweigh the benefit, but I'm optimistic.

The core of the idea is to determine all the plausible catalyst
placements before branching on them and then reuse this information as
much as possible in the subsequent branches. This would also give us
an opportunity to collect catalyst placements that "cause the same
perturbation", if we can come up with a good idea for what this means.

For each catalyst we will store a mask specifying whether the
placement is known to be valid and a mask for whether the placement
known to be not valid. Using these we can also remember that the
validity of a placement is unknown. "Valid" at least means that the
catalyst doesn't immediately explode, but there might be other useful
conditions.

We can step a configuration forward until a constraint is broken.
Let's call the cell and generation where this happens a "problem". The
most important constraints are when a catalyst is destroyed or a
filter is not matched.

At each node of the search tree:
1. Rollout the state and determine the problem with the current
   configuration, if any. If there is no problem, report the solution!
2. If there are no remaining interaction points in the current
   generation, advance.
2. Rollout (again) and collect the reaction envelope intersected
   with the lightcone shrinking towards the problem point. 
3. Determine which catalyst placements in that envelope are valid,
   using the cached validity masks where possible.
4. Collect the placements that cause the same perturbation of the
   active region and heuristically choose a good representative for
   each perturbation. (And remember the other placements for later.)
5. Branch on these perturbations. For each:
    1. Rollout (max 32 steps) and determine which interaction
       points' validity cannot be affected by this one. Reset the
       validity masks for all other placements.
    2. Recur!

When reporting a solution, try and complete the catalysts using
`cadical` (copy from Silk). If that fails, try all the other
combinations of catalyst choices

This is a bit rough and there are still some details to be worked out,
but I am hoping this will cut down massively on the number of
placements that have to be tested. As time advances, there are growing
lightcones from earlier catalyst placements and a shrinking lightcone
to the problem spot, and we only have to bother re-testing placements
that lie in the intersection.

For this to work well, I think we will need catalysts to specify more
information: the shapes of the active region at the spot that it
successfully interacts with a catalyst. For an eater, one of these
shapes might be

```
x = 7, y = 6, rule = LifeHistory
5.E$5.2E$2.2A$.A.A$.A$2A!
```

but there would also be others. These shapes might be larger for
catalysts with very specific reactions. And for transparent catalysts,
you could skip this check.

We could also allow a catalyst to specify an intermediate state that
it has to pass through to check it's recovering correctly. This would
have the biggest impact on traffic-stops and other things with
extremely long recovery times.

## Locked State

Rather than placing complete catalysts, we should just be placing the 

## TODO

Avoid catalyst placements that immediately interact with other catalysts

Properly collate perturbations done by different catalysts
    
Allow partial catalysts/weldling, by having a notion of board that fixes some of the cells to a certain state, can try welding at the end.

Copy over the stable completion code from silk

Do CountNeighbourhoodInteraction more intelligently rather than doing
the entire sum first.

Allow a single catalyst to have multiple approaches/reaction times

Automatically determine catalyst properties from some input soups

Option to report fizzles only (a la silk)

