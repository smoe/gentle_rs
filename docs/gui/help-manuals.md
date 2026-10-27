## Help manuals

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

The `Help` menu now includes:

The ordinary `Help` button in a window opens that window's manual topic. It is
separate from `Agent help`: reading the manual does not capture a screenshot,
send a model request or change the project. Sequence and specialist windows
retain the same documentation target in hosted and detached presentation.

The GUI manual is a book of [topic files](../gui.md#topics), also browsable on
GitHub. Its index is not a replacement for the full reference. Installed
resources are located relative to the application as well as the working
directory; a missing/incompatible external book falls back to all embedded
chapters together. Screenshots and tutorial artifacts still require installed
resources. Older `gui.md#section` links retain web entry points and resolve
directly to their topic in the application.

The CLI manual uses the same topic tree, section reader, complete-book search
and copy behavior. Related GUI/CLI topics link to each other near their heading;
Back and Forward retain navigation across both books. A topic without a direct
counterpart says so explicitly. These links aid navigation, not a claim that
every command has a graphical control. Old `cli.md#section` links also resolve
to their new chapters.

- `GUI Manual`: opens `docs/gui.md` one section at a time. Use **Section**,
  **Previous section** and **Next section** to navigate. Oversized sections
  have numbered continuations, split between complete Markdown blocks or
  top-level bullet items; code, tables, numbered lists and nested items stay together. Scrolling
  lays out only the selected section, not the entire manual. The section index
  and display text are cached until the source changes or **Reload** is used.
  **Find** searches the complete manual and selects the section containing
  the chosen match; the source line/snippet locates it within that section.
  **Selectable Text** shows the selected section, while **Copy Manual** copies
  the complete source. Other Help tabs and tutorials retain their own navigation.
  The GUI manual also offers a collapsible topic tree, clickable ancestor links
  and **Back** / **Forward** navigation. On narrow windows, expand **Contents**
  above the page; on wider windows the tree is alongside it. Selecting a heading
  link opens its section at the top, including same-page `#section` links and
  links from other help documents into `gui.md#section`. Missing GUI headings
  produce an explicit message. Manual navigation never invokes Agent Help.
- `CLI Manual`: opens `docs/cli.md` in an in-app markdown viewer
- `Agent Interface`: opens `docs/agent_interface.md` in an in-app markdown viewer
- `Reviewer Quickstart`: opens `docs/reviewer_preview.md` (internal preview checklist and known limitations)
- `Shell Commands`: generated command reference from `docs/glossary.json`
  (usage + task summary per command)
- `Tutorials`: opens tutorial markdown docs in the same help window
  (`Topic` selector starts the second help-header row, ahead of search)
  - curated order now comes from `docs/tutorial/catalog.json` when available
  - choose a topic group first, then its tutorial; catalog numbering and order
    are preserved, including entries discovered without catalog metadata
  - menus fit within the current window, with wrapped titles and scrollbars
    for long lists. Scroll with the wheel/trackpad, drag the scrollbar, or hold
    the pointer near the top/bottom of a list to scroll toward that end
  - falls back to recursive markdown discovery under `docs/tutorial/**` if the
    catalog is unavailable
  - discovery reads titles without preparing other tutorials' images. The
    document heading still takes precedence over the catalog title, with the
    same catalog/filename fallbacks. Only the selected tutorial prepares its
    images; the initial selected page remains preloaded with the manuals.
    Opening Help or the tutorial-project menu does not render the whole catalog.
    Images use the existing shared SVG font configuration and PNG cache, with
    the same original-image fallback if conversion fails; no network is needed.
  - when the selected catalog entry declares tutorial relationships, a compact
    row below the topic selector shows **Prerequisite** and **Continue with**
    buttons. Labels use the stable decimal tutorial id, hover text shows the
    full title and summary, and selecting a pointer reuses the same Help window.
    Each direction is limited to three explicit catalog ids; absent or fallback
    catalog metadata produces no guessed navigation.
  - inline and reference-style relative Markdown links also open the linked
    chapter in this Help window, resolved from the displayed document rather
    than the application's working directory. Existing catalog entries and
    manual tabs are reused; uncatalogued Markdown guides can be opened too.
    Missing targets produce a visible Help message instead of a silent browser
    request. External web links keep their normal browser handling.
    The following limitation applies outside the GUI and CLI manual books:
    A `chapter.md#section` link opens the page and discloses that automatic
    section scrolling is not supported; use Find to locate the section.
- on macOS, app menu `GENtle -> GENtle Help...` opens the same help window
- help now opens in its own native window (separate viewport), not as an overlay in the project window
- re-selecting an already open help page, including via `GENtle -> GENtle Windows`, reuses the existing help viewport instead of reloading the manuals from disk
- Shell command reference includes an `Interface` selector:
  - `All`
  - `GUI shell`
  - `CLI shell`
  - `CLI direct`
  - `JS`
  - `Lua`

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Open windows and focus](open-windows-and-focus.md) | [Next: Command Palette and History/Jobs Panels](command-palette-and-historyjobs-panels.md)
<!-- manual-page-nav:end -->
