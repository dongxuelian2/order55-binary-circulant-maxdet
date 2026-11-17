"""Replay comparison preserves contents across native newline conventions."""
from pathlib import Path
import runpy

api = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'scripts/verify_order55_global.py'))


def test_replay_newlines_and_corruption(tmp_path):
    saved, replay = tmp_path / 'saved', tmp_path / 'replay'
    saved.write_bytes(b'{"norm":17}\r\n{"norm":23}\r\n')
    replay.write_bytes(b'{"norm":17}\n{"norm":23}\n')
    assert api['same_replay'](saved, replay)
    replay.write_bytes(b'{"norm":17}\n{"norm":29}\n')
    assert not api['same_replay'](saved, replay)
    replay.write_bytes(b'{"norm":17}\n')
    assert not api['same_replay'](saved, replay)
