/* $Id$ */

/***
  This file is part of polypaudio.
 
  polypaudio is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published
  by the Free Software Foundation; either version 2 of the License,
  or (at your option) any later version.
 
  polypaudio is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with polypaudio; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA.
***/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "autoload.h"
#include "module.h"
#include "xmalloc.h"
#include "memchunk.h"
#include "sound-file.h"
#include "log.h"
#include "scache.h"

static void entry_free(struct pa_autoload_entry *e) {
    assert(e);
    pa_xfree(e->name);
    pa_xfree(e->module);
    pa_xfree(e->argument);
    pa_xfree(e->filename);
    pa_xfree(e);
}

static struct pa_autoload_entry* entry_new(struct pa_core *c, const char *name) {
    struct pa_autoload_entry *e = NULL;
    assert(c && name);
    
    if (c->autoload_hashmap && (e = pa_hashmap_get(c->autoload_hashmap, name)))
        return NULL;
    
    e = pa_xmalloc(sizeof(struct pa_autoload_entry));
    e->name = pa_xstrdup(name);
    e->module = e->argument = e->filename = NULL;
    e->in_action = 0;
    
    if (!c->autoload_hashmap)
        c->autoload_hashmap = pa_hashmap_new(pa_idxset_string_hash_func, pa_idxset_string_compare_func);
    assert(c->autoload_hashmap);
    
    pa_hashmap_put(c->autoload_hashmap, e->name, e);

    return e;
}

int pa_autoload_add_module(struct pa_core *c, const char*name, enum pa_namereg_type type, const char*module, const char *argument) {
    struct pa_autoload_entry *e = NULL;
    assert(c && name && module);

    if (!(e = entry_new(c, name)))
        return -1;
        
    e->module = pa_xstrdup(module);
    e->argument = pa_xstrdup(argument);
    e->type = type;
    return 0;
}

int pa_autoload_add_sample(struct pa_core *c, const char*name, enum pa_namereg_type type, const char*filename) {
    struct pa_autoload_entry *e = NULL;
    assert(c && name && filename);

    if (!(e = entry_new(c, name)))
        return -1;
        
    e->filename = pa_xstrdup(filename);
    e->type = PA_NAMEREG_SAMPLE;
    return 0;
}

int pa_autoload_remove(struct pa_core *c, const char*name, enum pa_namereg_type type) {
    struct pa_autoload_entry *e;
    assert(c && name && type);

    if (!c->autoload_hashmap || !(e = pa_hashmap_get(c->autoload_hashmap, name)))
        return -1;

    pa_hashmap_remove(c->autoload_hashmap, e->name);
    entry_free(e);
    return 0;
}

void pa_autoload_request(struct pa_core *c, const char *name, enum pa_namereg_type type) {
    struct pa_autoload_entry *e;
    struct pa_module *m;
    assert(c && name);

    if (!c->autoload_hashmap || !(e = pa_hashmap_get(c->autoload_hashmap, name)) || (e->type != type))
        return;

    if (e->in_action)
        return;

    e->in_action = 1;

    if (type == PA_NAMEREG_SINK || type == PA_NAMEREG_SOURCE) {
        if ((m = pa_module_load(c, e->module, e->argument)))
            m->auto_unload = 1;
            
    } else {
        struct pa_sample_spec ss;
        struct pa_memchunk chunk;
        assert(type == PA_NAMEREG_SAMPLE);

        if (pa_sound_file_load(e->filename, &ss, &chunk, c->memblock_stat) < 0) 
            pa_log(__FILE__": failed to load sound file '%s' for autoload entry '%s'.\n", e->filename, e->name);
        else {
            pa_scache_add_item(c, e->name, &ss, &chunk, NULL, 1);
            pa_memblock_unref(chunk.memblock);
        }
    }

    e->in_action = 0;
}

static void free_func(void *p, void *userdata) {
    struct pa_autoload_entry *e = p;
    entry_free(e);
}

void pa_autoload_free(struct pa_core *c) {
    if (!c->autoload_hashmap)
        return;

    pa_hashmap_free(c->autoload_hashmap, free_func, NULL);
}
