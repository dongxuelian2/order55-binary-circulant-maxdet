# Proof dependency graph

```text
M15 transcription + μ₃ ring arithmetic
        │
        ├── exact benchmark norm D_record
        ├── inner-product catalogue ──┐
        │                             ├── trace stability: Q ≤ 168
        │                             └── norm-mod-9 color congruence
        │                                  │
        │                                  ├── F7 exclusion
        │                                  ├── code independence B₃(15,10)=12
        │                                  └── four color partitions
        │                                       │
        │                    ┌──────────────────┴──────────────────┐
        │                    │                                     │
        │          exact support/component certificates       Schur/spectral caps
        │                    │                                     │
        │                    └──────────────┬──────────────────────┘
        │                                   │
        │                    Q<144 exact aggregate (replay passed)
        │                                   │
        │                         Q=144 boundary ── Q=150 boundary
        │                                   │             │
        │                                   └──────┬──────┘
        │                                          │
        └── Q=153,159,162,168 trace envelope; Q≥171 below record
                                                   │
                              3^14 divisibility + Eisenstein norm sieve
                                                   │
                                current rigorous integer upper (not equality)
```

The dashed conceptual endpoint is deliberate: no branch currently proves
that the upper bound is attained or that the genuine Q=153,159,162,168 trace
envelope is realizable.
