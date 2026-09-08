# GLOSSARY — ffpopt

> Canonical definitions for domain terms used in this repo. The aim is
> *consistency*: one definition per term, one place to update it when the
> meaning shifts.

## What goes in

- Domain terms with non-obvious meaning (a reader who has not worked on
  ffpopt for six months should reach for this file).
- Terms used in three or more docs or modules with a specific repo meaning
  that diverges from the generic English meaning.
- Cross-cutting concepts that appear in `dev/claude-skills/` skill files,
  README, and code identifiers — anchor them once here, point everywhere else.

## What does NOT go in

- Private internal names (variables, classes) — those live in the code.
- Generic programming concepts (database, queue, retry) unless this repo uses
  them in a non-standard way.
- Per-team jargon that belongs in a chat channel description, not here.

## How to add an entry

Use the shape below. Lead with the name. One-sentence definition. Add detail
only when the one-sentence form is genuinely ambiguous or invites misuse.
Always cite an authoritative source — a file path plus section.

```markdown
### TermName

**Definition.** One-sentence definition.

**Detail.** Optional longer explanation, including when the term applies and
common confusions with similar terms.

**Authoritative source.** `<path>:<section or line range>`
```

## Platform primitives

<entries for repo-wide structural concepts — orchestration objects, core
types, lifecycle states>

## Domain concepts

<entries for scientific, business, or product-domain terms — the vocabulary
the team uses when discussing the work, not the code>

## Operational terms

<entries for deploy targets, environments, runtime concepts that span
multiple docs>

## External integrations

<entries for third-party services, APIs, providers as named in this repo
specifically — e.g., how this repo refers to "Atlas" vs the generic term>
