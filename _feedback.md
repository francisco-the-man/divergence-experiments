# Divergence — feedback queue

_Source of truth: the SQLite ticket DB on the Modal Volume. Flushed every 5 min by `flush_feedback_tick`. To mark an item addressed, call `mark_feedback_addressed(<id>, "what you did")` via the MCP tools, or update the DB directly._

**0 open** · 1 recently addressed (last 10)

## Open

_No open feedback. Inbox zero._

---

## Recently addressed

### `fb_02de2be4` — *meta* — addressed 2026-05-18T17:55:10+00:00

> TEST: round-trip check from claude-code session. delete me whenever.

**Resolution:** round-trip verified; pipeline operational, see commit f679180+
